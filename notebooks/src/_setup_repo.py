# --- Find the MEGAN_SC folder and, in Google Colab, connect Google Drive (the same block opens every notebook) ---------
import os, sys, glob
try:
    from google.colab import drive as _drive
    IN_COLAB = 'MEGAN_SC_RUNNER' not in os.environ          # tools/run_notebook.py sets it when it runs a notebook outside Colab
except ImportError:
    IN_COLAB = False
if IN_COLAB:
    # A failed mount must stop the notebook: otherwise it would read and write a local /content/drive folder that is not
    # Google Drive.
    if os.path.isdir('/content/drive') and not os.path.ismount('/content/drive') and os.listdir('/content/drive'):
        raise RuntimeError('/content/drive holds local files written while Google Drive was not mounted, so Drive cannot be '
                           'mounted there. Use Runtime > Disconnect and delete runtime, then run the notebook again.')
    _drive.mount('/content/drive')


def find_repository():
    """The MEGAN_SC folder: $MEGAN_SC_REPO when set; in Colab, MyDrive/MEGAN_SC (or a MEGAN_SC folder one level down in
    MyDrive, or /content/MEGAN_SC for a clone outside Drive); elsewhere, the current folder or its parent."""
    cands = [os.environ.get('MEGAN_SC_REPO', '')]
    if IN_COLAB:
        cands += ['/content/drive/MyDrive/MEGAN_SC', '/content/drive/MyDrive/megan_sc', '/content/drive/MyDrive/megan-sc']
        cands += sorted(glob.glob('/content/drive/MyDrive/*/MEGAN_SC')) + ['/content/MEGAN_SC']
    cands += [os.getcwd(), os.path.dirname(os.getcwd())]
    for c in cands:
        if c and os.path.isfile(os.path.join(c, 'src', 'megan_sc', 'fixes.py')):
            return os.path.abspath(c)
    raise SystemExit('The MEGAN_SC folder was not found. In Colab, copy the repository to Google Drive as MyDrive/MEGAN_SC '
                     '(tutorial, chapter 3), or run  os.environ["MEGAN_SC_REPO"] = "/path/to/MEGAN_SC"  before this cell.')


REPO_DIR = find_repository()
if os.path.join(REPO_DIR, 'src') not in sys.path:
    sys.path.insert(0, os.path.join(REPO_DIR, 'src'))
print('MEGAN_SC folder:', REPO_DIR, '| Colab:', IN_COLAB)
