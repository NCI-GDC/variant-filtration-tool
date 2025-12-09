from unittest import TestCase
from unittest.mock import MagicMock, Mock, call, mock_open, patch

from gdc_filtration_tools.tools.format_strelka_vcf import (
    add_filter,
    adjust_INDEL,
    adjust_record,
    adjust_SNV,
    convert_gt_spec,
    ensure_gt,
    format_strelka_vcf,
    parse_field,
    parse_info,
    qsi_filter,
)


class TestFormatStrelka(TestCase):
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.print")
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.add_filter")
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.ensure_gt")
    @patch(
        "gdc_filtration_tools.tools.format_strelka_vcf.adjust_record",
        return_value="new_row",
    )
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.tabix_index")
    def test_format_strelka_vcf(
        self, tabix_index, adjust_record, ensure_gt, add_filter, print_fn
    ):
        vcf = MagicMock()
        # fd = {"FORMAT": 'header_format'}
        vcf.header.__getitem__.side_effect = ["header_format", "header_filter"]
        vcf.iter_header_lines = Mock(return_value=["header_line"])
        vcf.iter_rows = Mock(return_value=["row"])
        open_fn = mock_open()

        with (
            patch("gdc_filtration_tools.tools.format_strelka_vcf.open", open_fn),
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.VcfReader",
                return_value=vcf,
            ) as vcfreader,
        ):
            format_strelka_vcf("input.vcf", "out.vcf.gz")

            vcfreader.assert_called_once_with("input.vcf")
            ensure_gt.assert_called_once_with("header_format")
            add_filter.assert_called_once_with("header_filter")
            open_fn.assert_called_once_with("out.vcf", "wt")
            handle = open_fn()
            print_call_list = [
                call("header_line", file=handle),
                call("new_row", file=handle),
            ]
            assert print_fn.call_args_list == print_call_list
            adjust_record.assert_called_once_with("row")
            tabix_index.assert_called_once_with("out.vcf", preset="vcf")

    def test_ensure_gt(self):
        fs = {"FOO": "FOO"}
        expected = fs.copy()
        expected["GT"] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        result = ensure_gt(fs)
        assert result == expected

    def test_add_filter(self):
        fs = {"FOO": "FOO"}
        expected = fs.copy()
        expected["LowQSI"] = (
            '##FILTER=<ID=LowQSI,Description="QSI value is at or below 10">'
        )
        result = add_filter(fs)
        assert result == expected

    def test_adjust_SNV(self):
        row = Mock()
        row.NORMAL = "3:5"
        row.TUMOR = "4:6"
        row.INFO = "row_info"
        row.FORMAT = "row_format"

        with (
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
                return_value={"NT": "nt value"},
            ) as pi,
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.convert_gt_spec",
                return_value="0/0",
            ) as cgs,
        ):
            adjust_SNV(row)
            pi.assert_called_once_with(row.INFO)
            cgs.assert_called_once_with("nt value")
            row.replace.assert_called_once_with(
                NORMAL="0/0:3:5", TUMOR="0/1:4:6", FORMAT="GT:row_format"
            )

    def test_adjust_INDEL(self):
        row = Mock()
        row.NORMAL = "3:5"
        row.TUMOR = "4:6"
        row.INFO = "row_info"
        row.FORMAT = "row_format"
        row.FILTER = "PASS"

        with (
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
                return_value={"NT": "germline", "SGT": "foo->somatic", "QSI": 5},
            ) as pi,
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.convert_gt_spec",
                side_effect=["0/0", "0/1"],
            ) as cgs,
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.qsi_filter",
                return_value=True,
            ),
        ):
            adjust_INDEL(row)
            pi.assert_called_once_with(row.INFO)
            assert cgs.call_args_list == [call("germline"), call("somatic")]
            row.replace.assert_called_once_with(
                NORMAL="0/0:3:5",
                TUMOR="0/1:4:6",
                FORMAT="GT:row_format",
                FILTER="LowQSI",
            )

    @patch(
        "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
        return_value={"QSI": "5"},
    )
    def test_qsi_filter(self, pi):
        row = Mock()
        row.INFO = "row_info"
        result = qsi_filter(row)
        pi.assert_called_once_with("row_info")
        assert result

    def test_convert_gt_spec(self):
        assert convert_gt_spec("ref") == "0/0"
        assert convert_gt_spec("het") == "0/1"
        assert convert_gt_spec("hom") == "1/1"
        assert convert_gt_spec("conflict") == "./."

    @patch("gdc_filtration_tools.tools.format_strelka_vcf.adjust_INDEL")
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.adjust_SNV")
    def test_adjust_record(self, adjust_snv, adjust_indel):
        indel_info_key_set = {
            "IC",
            "IHP",
            "QSI",
            "OVERLAP",
            "QSI_NT",
            "RC",
            "RU",
            "TQSI",
            "TQSI_NT",
        }
        snv_info_key_set = {
            "ACGTNacgtnMINUS",
            "ACGTNacgtnPLUS",
            "DP",
            "QSS",
            "QSS_NT",
            "ReadPosRankSum",
            "SNVSB",
            "TQSS",
            "TQSS_NT",
        }
        common_key_set = {
            "MQ",
            "MQ0",
            "NT",
            "SGT",
            "SOMATIC",
            "SomaticEVS",
        }
        indel_set = indel_info_key_set | common_key_set
        snv_set = snv_info_key_set | common_key_set
        row = Mock()
        indel_keydict = {k: v for k, v in zip(indel_set, range(len(indel_set)))}
        snv_keydict = {k: v for k, v in zip(snv_set, range(len(snv_set)))}
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=indel_keydict,
        ):
            adjust_record(row)
            adjust_indel.assert_called_once_with(row)

        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=snv_keydict,
        ):
            adjust_record(row)
            adjust_snv.assert_called_once_with(row)

    def test_adjust_record_unknown(self):
        key_set = {"NOT" "EXPECTED"}
        row = Mock()
        info_keydict = {k: v for k, v in zip(key_set, range(len(key_set)))}
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=info_keydict,
        ):
            self.assertRaises(ValueError, adjust_record, row)

    @patch(
        "gdc_filtration_tools.tools.format_strelka_vcf.parse_field",
        side_effect=[("one", 1), ("two", 2)],
    )
    def test_parse_info(self, pf):
        result = parse_info("one;two")
        assert pf.call_args_list == [call("one"), call("two")]
        assert result == {"one": 1, "two": 2}

    def test_parse_field(self):
        result = parse_field("key=value")
        assert result == ("key", "value")

        result = parse_field("key")
        assert result == ("key", True)
