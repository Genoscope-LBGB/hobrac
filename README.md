# hobrac

hobrac is a rdbioseq project

A new version of the environment module associated with this project is deployed upon every commit on master, except if the commit message contains [skip ci]. The default version increment is patch. To ship a new minor or major version of the module, update moduleversion in .gitlab-ci.yml

It is recommanded that you setup tests for hobrac in as many `test_*.py` scripts as you want inside the `tests` folder: they will be executed by `pytest` everytime you push to GitLab, and the deployment stage will not be triggered unless they are successful

This lib can be installed using `pip install`. Once installed, the scripts in `bin/` will be available in your $PATH
