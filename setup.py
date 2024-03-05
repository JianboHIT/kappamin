from setuptools import setup

setup(
    name='kappamin',
    version='0.2.1',
    url='https://github.com/JianboHIT/kappamin',
    author='Jianbo ZHU',
    description='A python3 code for calculations of the minimum limit to thermal conductivity',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    keywords=['thermal-conductivity', 'limitation', 'condensed-matter-physics', 'phonon-transport'],
    license='Apache-2.0 license',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
    ],
    py_modules=['kappamin'],
    install_requires=[
        'numpy',
        'scipy'
    ],
)
