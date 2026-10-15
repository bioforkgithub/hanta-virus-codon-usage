import sys
print("python", sys.version.split()[0])
try:
    import Bio
    print("biopython", Bio.__version__)
except Exception as e:
    print("biopython MISSING", e)
try:
    import scipy
    print("scipy", scipy.__version__)
except Exception as e:
    print("scipy MISSING", type(e).__name__)
try:
    import urllib.request
    with urllib.request.urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        "?db=nuccore&term=%22Rockport+orthohantavirus%22%5BOrganism%5D&retmax=0",
        timeout=30,
    ) as r:
        print("NET OK", r.read().decode()[:400])
except Exception as e:
    print("NET FAIL", type(e).__name__, e)
