import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="xmascat",
    version="0.1.0",
    install_requires=[
        "pandas",
	"xarray",
	"astropy",
	"xarray_dataclasses", 
	"python-casacore", 
	"netCDF4",
	"joblib",
    ],
#    entry_points={
#        'console_scripts': [
#            'xmascat=xmascat:main',
#        ],
#    },
    author="ShinjiFujita",
    author_email="fujita.shinji.0728@gmail.com",
    description="Xarray to Measurementset(v2) conversion module for ASte and other Common Astronomical Telescopes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShinjiFujita/xmascat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
