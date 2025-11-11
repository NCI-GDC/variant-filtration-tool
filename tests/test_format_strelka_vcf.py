from unittest import TestCase
from unittest.mock import MagicMock, Mock, call, patch

from gdc_filtration_tools.tools.format_strelka_vcf import (
    adjust_INDEL,
    adjust_record,
    adjust_SNV,
    convert_gt_spec,
    ensure_gt,
    format_strelka_vcf,
    get_indel_or_snp_fn,
    parse_field,
    parse_info,
)


class TestFormatStrelka(TestCase):
    @patch("__main__.open")
    @patch("__main__.print")
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.VcfReader")
    @patch("gdc_filtration_tools.tools.format_strelka_vcf.ensure_gt")
    def test_format_strelka_vcf(self, ensure_gt, vcfreader, print_fn, open_fn):
        format_strelka_vcf("input.vcf", "output.vcf")

        vcfreader.assert_called_once_with("input.vcf")
        ensure_gt.assert_called_once()

    def test_ensure_gt(self):
        fs = {"FOO": "FOO"}
        expected = fs.copy()
        expected["GT"] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        result = ensure_gt(fs)
        assert result == expected

    def test_adjust_record(self):
        adjust_fn = Mock()
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.get_indel_or_snp_fn",
            return_value=adjust_fn,
        ) as iors:
            row = ("one", "two")
            adjust_record(row)
            iors.assert_called_once_with(row)
            adjust_fn.assert_called_once_with(row)

    def test_adjust_SNV(self):
        row = Mock()
        row.NORMAL = "bad:3:5"
        row.TUMOR = "4:6:bad"
        row.INFO = "row_info"

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
            row._replace.assert_called_once_with(NORMAL="0/0:3:5", TUMOR="0/1:4:6")

    def test_adjust_INDEL(self):
        row = Mock()
        row.NORMAL = "bad:3:5"
        row.TUMOR = "4:6:bad"
        row.INFO = "row_info"

        with (
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
                return_value={"NT": "germline", "SGT": "foo->somatic"},
            ) as pi,
            patch(
                "gdc_filtration_tools.tools.format_strelka_vcf.convert_gt_spec",
                side_effect=["0/0", "0/1"],
            ) as cgs,
        ):
            adjust_INDEL(row)
            pi.assert_called_once_with(row.INFO)
            assert cgs.call_args_list == [call("germline"), call("somatic")]
            row._replace.assert_called_once_with(NORMAL="0/0:3:5", TUMOR="0/1:4:6")

    def test_convert_gt_spec(self):
        assert convert_gt_spec("ref") == "0/0"
        assert convert_gt_spec("het") == "0/1"
        assert convert_gt_spec("hom") == "1/1"
        assert convert_gt_spec("conflict") == "./."

    def test_get_indel_or_snp_fn_indel(self):
        indel_set = {
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
        common_key_set = {
            "MQ",
            "MQ0",
            "NT",
            "SGT",
            "SOMATIC",
            "SomaticEVS",
        }
        all_keys = indel_set | common_key_set
        row = Mock()
        info_keydict = {k: v for k, v in zip(all_keys, range(len(all_keys)))}
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=info_keydict,
        ) as pi:
            result = get_indel_or_snp_fn(row)
            assert result is adjust_INDEL

    def test_get_indel_or_snp_fn_snv(self):
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
        all_keys = snv_info_key_set | common_key_set
        row = Mock()
        info_keydict = {k: v for k, v in zip(all_keys, range(len(all_keys)))}
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=info_keydict,
        ) as pi:
            result = get_indel_or_snp_fn(row)
            assert result is adjust_SNV

    def test_get_indel_or_snp_fn_error(self):
        key_set = {"NOT" "EXPECTED"}
        row = Mock()
        info_keydict = {k: v for k, v in zip(key_set, range(len(key_set)))}
        with patch(
            "gdc_filtration_tools.tools.format_strelka_vcf.parse_info",
            return_value=info_keydict,
        ) as pi:
            self.assertRaises(ValueError, get_indel_or_snp_fn, row)

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
