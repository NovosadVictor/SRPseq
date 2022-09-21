import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='SRPseq',
    version='0.0',
    scripts=[
        'srpseq',
    ],
    author='HSE.Bioinformatics',
    author_email='victor.o.novosad@gmail.com',
    description='Splicing Regulation Prediction from RNA-seq data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/NovosadVictor/SRPseq',
    packages=setuptools.find_packages(),
    install_requires=[
        'scipy',
        'scikit-learn',
        'numpy',
        'pandas',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
    ],
)
