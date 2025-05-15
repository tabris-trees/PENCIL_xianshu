from setuptools import setup, find_packages

setup(
    name="pencil",
    version="2.0",
    description="Python for pencil postprocessing scripts",
    # packages=find_packages(),
    packages=[
        "pencil",
        "pencil.backpack",
        "pencil.calc",
        "pencil.diag",
        "pencil.export",
        "pencil.io",
        "pencil.ism_dyn",
        "pencil.math",
        "pencil.read",
        "pencil.sim",
        "pencil.tool_kit",
        "pencil.util",
        "pencil.visu",
    ],
)
