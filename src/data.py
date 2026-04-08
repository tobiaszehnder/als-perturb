import pickle
import urllib.request
import anndata as ad
from pathlib import Path

DATA_URL = (
    "https://s3.eu-west-2.amazonaws.com/helical-candidate-datasets/counts_combined_filtered_BA4_sALS_PN.h5ad"
    "?X-Amz-Algorithm=AWS4-HMAC-SHA256"
    "&X-Amz-Credential=AKIA4SYAMZBOLGEU2HNG%2F20260401%2Feu-west-2%2Fs3%2Faws4_request"
    "&X-Amz-Date=20260401T180147Z"
    "&X-Amz-Expires=604800"
    "&X-Amz-SignedHeaders=host"
    "&X-Amz-Signature=ed35cf3e999a7ea692fe3127a87a25d71f180dbded67e5522aee0d78a6cb4a7b"
)

DEFAULT_DATA_PATH = Path(__file__).parent.parent / "data" / "counts_combined_filtered_BA4_sALS_PN.h5ad"
DEFAULT_EMBEDDINGS_PATH = Path(__file__).parent.parent / "data" / "embeddings.pkl"


def load_embeddings(path: Path = DEFAULT_EMBEDDINGS_PATH) -> dict:
    path = Path(path)
    if path.exists():
        with open(path, "rb") as f:
            embeddings = pickle.load(f)
        print(f"Loaded {len(embeddings)} existing embeddings from {path}: {list(embeddings.keys())}")
        return embeddings
    print("No existing embeddings found, starting fresh.")
    return {}


def save_embeddings(embeddings: dict, path: Path = DEFAULT_EMBEDDINGS_PATH) -> None:
    path = Path(path)
    path.parent.mkdir(exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(embeddings, f)


def load_data(path: Path = DEFAULT_DATA_PATH) -> ad.AnnData:
    """Download (if needed) and load the ALS dataset."""
    path = Path(path)
    path.parent.mkdir(exist_ok=True)

    if not path.exists():
        print("Downloading dataset...")
        urllib.request.urlretrieve(DATA_URL, path)
        print("Download complete.")

    adata = ad.read_h5ad(path)
    print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata
