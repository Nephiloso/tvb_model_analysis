from setuptools import setup, find_packages

setup(
      name = 'Crit',
      version = 1.0.0,
      url ='',
      description = '',
      packages = find_packages
      install_requires = [
          Crit @ git +ssh://git@github.com/nanda2502/Crit.git]
)
