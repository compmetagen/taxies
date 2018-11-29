from setuptools import setup, find_packages

setup(
    name='taxies',
    version='0.1.dev',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        taxies=taxies.scripts.taxies:cli
    ''',
)