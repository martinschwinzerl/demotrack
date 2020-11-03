import setuptools
import sys

requirements = {"install": ["numpy", "scipy", "pyopencl"]}

if sys.version_info < (3, 7):
    requirements["install"].append("dataclasses")

version = open("demotrack/__init__.py").readline().split('"')[1]

setuptools.setup(
    name="demotrack",
    version=version,
    description="Simplified 6D demo tracking code",
    author="Martin Schwinzerl",
    author_email="martin.schwinzerl@cern.ch",
    url="https://github.com/martinschwinzerl/demotrack",
    packages=["demotrack", ],
    package_dir={"demotrack": "demotrack"},
    package_data = {
        "demotrack": [ "src/beam_elements.h", "src/track.h", "src/particle.h" ] },
    install_requires=requirements["install"],
)
