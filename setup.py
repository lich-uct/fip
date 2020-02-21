
#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='fip',
      version='0.0.1',
      description='Python library for profiling interrelations between qualitative features',
      url='https://github.com/lich-uct/fip',
      author='Ivan Cmelo',
      license='GPL-3.0',
      packages=find_packages(),
      test_suite='test',
      include_package_data=True
     )
