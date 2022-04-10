import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FastAAI",
    version="0.1.08",
    author="Kenji Gerhardt",
    author_email="kenji.gerhardt@gmail.com",
    description="A rapid estimator for amino acid identity between genomes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
	#py_modules=["fastaai"],
	include_package_data=True,
    python_requires='>=3',
	entry_points={
        "console_scripts": [
            "fastaai=fastaai.FastAAI:main",
        ]
    }
)

