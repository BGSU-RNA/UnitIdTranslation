from setuptools import setup

setup(
    name = "id_translate",
    version = "0.1.4",
    author = "Blake Sweeney",
    author_email = "bsweeney@bgsu.edu",
    description = ("Translate old to new style ids"),
    license = "BSD",
    py_modules = ["id_translate"],
    requires = ['biopython>=1.60'],
)
