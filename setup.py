import setuptools

setuptools.setup(
    name="mandoline",
    packages=["mandoline"],
    version="0.0.0",
    install_requires=["sortedcontainers", "colorama"],
    entry_points={
        "console_scripts": [
            "mandoline_build_counting_dag=mandoline:build_counting_dag",
            "mandoline_count=mandoline:count",
            "mandoline_decomposition=mandoline:decomposition",
            "mandoline_show_pattern_decomposition=mandoline:show_pattern_decomposition",
            "mandoline_enumerate=mandoline:enumerate_",
        ]
    },
)
