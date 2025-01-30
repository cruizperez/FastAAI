import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="FastAAI-release",
	version="0.0.01",
	author="Kenji Gerhardt",
	author_email="kenji.gerhardt@gmail.com",
	description="A rapid estimator for amino acid identity between genomes.",
	long_description=long_description,
	long_description_content_type="text/markdown",
	packages=setuptools.find_packages(),
	py_modules=["fastaai_api"],
	include_package_data=True,
	python_requires='>=3',
	install_requires=[
		'numpy==2.2.2',
		'pyrodigal==3.6.3',
		'pyhmmer==0.11.0',
	],
	entry_points={
		"console_scripts": [
			"fastaai=fastaai.fastaai:main",
		]
	}
)

