import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dimcheck",  # 模块名称
    version="1.0.2",  # 当前版本
    author="Tony Riddick",  # 作者
    author_email="tonyriddick1527@gmail.com",  # 作者邮箱
    description="A physics dimension checker based on sympy",  # 模块简介
    long_description=long_description,  # 模块详细介绍
    long_description_content_type="text/markdown",  # 模块详细介绍格式
    url="https://github.com/Tony1527/dimcheck",  # 模块github地址
    packages=setuptools.find_packages(exclude=[".history"],include=["dimcheck"]),  # 自动找到项目中导入的模块
    include_package_data=True,
    readme="README.md",
    exclude_package_data={"":[".gitignore",".history","__pycache__"]},
    package_data={"dimcheck":["*.json"]},
    # 模块相关的元数据
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    keywords=["physics", "dimension", "dimension analysis", "unit checker", "dimension checker", "sympy","量纲分析","量纲检查","单位检查","量纲","单位","sympy"],
    # 依赖模块
    install_requires=[
        'sympy',
        'numpy',
    ],
    python_requires='>=3.7'
)
