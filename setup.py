import setuptools

setuptools.setup(
    name="mandoline",
    packages=["mandoline"],
    version="0.1.0",
    install_requires=["sortedcontainers", "colorama"],
    entry_points={
        "console_scripts": [
            "mandoline=mandoline.__main__:main"
            # "mandoline_decompose=mandoline:build_counting_dag",
            # "mandoline_count=mandoline:count",
            # "mandoline_enumerate=mandoline:enumerate_",
        ]
    },
)
