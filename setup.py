import setuptools

setuptools.setup(
    name= 'meltInc',
    version= '0.1',
    description= '...',
    author= 'Thomas van Gerve',
    
    packages= setuptools.find_packages(
        exclude= ['examples']
        ),

    # package_dir= {'' : 'petroPy'},
    package_data= {'meltInc': ['static/*']}, 

    install_requires= [
    'pandas',
    'matplotlib',
    'numpy',
    'scipy',
    'sympy',
    'csaps'
    ]
)