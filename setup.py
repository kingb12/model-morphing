from setuptools import setup

setup(name='model_morphing',
      version='0.1',
      description='A python package for the model morphing workflow',
      author='Brendan King <bking@systemsbiology.org',
      author_email='bking@systemsbiology.org',
      license='MIT',
      packages=['lib', 'test'],
      install_requires=[
          'cobra',
      ],
      zip_safe=False)