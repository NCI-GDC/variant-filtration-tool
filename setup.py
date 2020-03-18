import setuptools

setuptools.setup(
    name = "gdc-filtration-tools",
    author = "Kyle Hernandez",
    author_email = "kmhernan@uchicago.edu",
    version = 1.0,
    description = "Variant filtration utilities for GDC workflows.",
    url = "https://github.com/NCI-GDC/variant-filtration-tool",
    license = "Apache 2.0",
    packages = setuptools.find_packages(),
    entry_points = {
        'console_scripts': [
            'gdc-filtration-tools = gdc_filtration_tools.__main__:main'
        ]
    },
    install_requires = [
      'pysam',
      'defopt'
    ]
)
