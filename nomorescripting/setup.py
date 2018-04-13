from setuptools import setup




setup(name='NoMoreScripting',
      version='0.1',
      description='Set of python scripts to make my bioinformatics life easier',
      url='https://github.com/Jome0169/NoMoreScripting',
      author='Pablo Mendieta',
      author_email='pablo.mendieta75@gmail.com',
      license='MIT',
      packages=['matplotlib', 'numpy', 'BioPython', 'statistics' ],
      install_requires=[
          'python3',
      ],
      zip_safe=False)



