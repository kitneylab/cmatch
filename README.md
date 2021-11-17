# cMatch 

## Prerequisite

You need Python 3.9 and pip the package installer for Python [pip](https://pip.pypa.io/en/stable/)


## Install

Install a virtual environment if you want.

```
    $ python3 -m venv --prompt cmatch venv
    $ source ./venv/bin/activate
```

Install teh dependencies

```
    $ pip install -r requirements.txt
```

## cMatch command line tool

run cmatch.py

```
    $ ./cmatch.py --help 
```

## Simple Example 

Run cMatch on a simple example, matching the violacein template and two synthetic violacein constructs

```
    $ cd simple_example
    $ ../cmatch.py template.json vio-B0030-B0030-B0030-B0030-B0030.seq vio-B0030-B0031-B0032-B0033-B0064.seq
```

the tool will output a JSON


## Run the tests

### Run the tests for algorithm CM_0

```
    $ python Testing_Algorithm_CM_0.py
```

### Run the tests for algorithm CM_1

```
    $ python Testing_Algorithm_CM_1.py
```

### Run the tests for algorithm CM_2 First Example

```
    $ python test_cm2_vio_easy_1vs1_th75.py
    $ python test_cm2_vio_easy_1vs1_th99.py
    $ python test_cm2_vio_hard_1vs1_th75.py
    $ python test_cm2_vio_hard_1vs1_th99.py
```

### Run the tests for algorithm CM_2 Second Example

```
    $ python test_cm2_vio_easy_1vsAll_th75.py
    $ python test_cm2_vio_easy_1vsAll1_th99.py
    $ python test_cm2_vio_hard_1vsAll1_th75.py
    $ python test_cm2_vio_hard_1vsAll1_th99.py
```

### Run the tests for algorithm CM_2 Real life example Lycopene Operon 

```
    $ python test_cm2_lycopene_sanger_10.py
    $ python test_cm2_lycopene_sanger_100.py

```

### Run the tests for algorithms CM_1, CM_2 Violacein-0000 cat x2 and cat x3

```
    $ python cm1_cat2.py
    $ python cm1_cat2min.py
```

```
    $ python cm1_cat3.py
    $ python cm1_cat3min.py
```

```
    $ python cm2_cat2.py
    $ python cm2_cat2min.py
```

```
    $ python cm2_cat3.py
    $ python cm2_cat3min.py
```
