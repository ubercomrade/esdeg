# Tool name

## Introduction
Text

## Python requirements

  * numpy
  * scipy
  * numba
  * matplotlib

  You can use pip to install dependences `pip install numpy, scipy, numba, matplotlib`,
  Or you can use conda `conda install numpy, scipy, numba, matplotlib`

## Installation

```  
git clone https://github.com/XXX/XXX.git  
cd pipeline/  
pip install -e .  
```

## Usage
The command `enrest.py -h` return:

```
test
```

### Example run
```
None
```

### Required arguments description

**First positional argument**:
```
None
```
PASS



**Second positional argument**:
```
N                     promoters of organism (hg38, mm10, tair10, b73)
```
Value of N can be _hg38_ or _mm10_. Depend on organism usage in research

**Third positional argument**:
```
None
```
PASS

**Fourth positional argument**:
```
output                output dir
```
Path of directory to write results. If dir is not exists it'll be created.


### Optional arguments description

```
None
```
PASS

```
None
```
PASS

## License

Copyright (c) 2021 Anton Tsukanov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.