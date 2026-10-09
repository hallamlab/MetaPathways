import importlib.util
from pathlib import Path
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('compile_recipe', ROOT / 'conda_recipe/compile_recipe.py')
recipe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recipe)


class RecipeTests(unittest.TestCase):
    def test_source_environment_pip_does_not_become_runtime_dependency(self):
        dependencies = yaml.safe_load((ROOT / 'docker/conda_base.yml').read_text())['dependencies']
        self.assertIn('pip', dependencies)
        runtime = recipe.runtime_dependencies(dependencies)
        self.assertNotIn('pip', runtime)
        self.assertIn('nextflow>=25.10,<27', runtime)
        self.assertIn('urllib3>=2.8.0,<3', runtime)
        self.assertIn('setuptools>=83,<85', runtime)
        template = yaml.safe_load((ROOT / 'conda_recipe/meta_template.yaml').read_text().replace('<HELPER_SOURCES>', '').replace('<ENTRY>', '').replace('<REQUIREMENTS>', ''))
        self.assertIn('pip', template['requirements']['host'])

    def test_only_pip_is_filtered_and_vcs_entries_are_rejected(self):
        self.assertEqual(recipe.runtime_dependencies(['pip>=26', 'conda-forge::pip=26', 'pip-tools', 'python=3.11']), ['pip-tools', 'python=3.11'])
        with self.assertRaisesRegex(ValueError, 'explicit package strings'):
            recipe.runtime_dependencies([{'pip': ['some-vcs-source']}])
