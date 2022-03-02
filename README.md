# Evolution of evolvability in changing environments

(working title!)


## Dependencies

Python dependencies for running job scripts + data aggregation/management scripts:

```bash
python3 -m venv pyenv
source pyenv/bin/activate
pip3 install -r requirements.txt
```

## Building Avida

This repository contains the version of Avida that we used to run our experiments.

From a fresh clone of this repository, first you'll need to init and update the git submodules:

```
git submodule init
git submodule update
```

Next, you should be able to run the `build_avida` script:

```
cd avida
./build_avida
```