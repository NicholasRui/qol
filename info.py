import qol
import subprocess

# main qol directory
qol_path = qol.__file__.replace('__init__.py', '')

# attempt to get qol commit hash
try:
    qol_commit_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL).decode('utf-8').strip()
except:
    qol_commit_hash = 'unknown'
