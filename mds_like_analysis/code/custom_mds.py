import numpy as np
from scipy.spatial.distance import pdist, squareform, cdist


class CustomMDS:
    """
    Implements Classical MDS with approximate out-of-sample extension.

    - fit(X): computes the classical MDS embedding for the training data X
    - transform(X_new): embeds new points approximating the out-of-sample extension
    - fit_transform(X): combines fit and transform in a single step
    """

    def __init__(self, n_components=2):
        self.n_components = n_components
        self.embedding_ = None
        self.X_fit_ = None
        self.eigvecs_ = None
        self.L_sqrt_ = None

    def fit(self, X):
        X = np.asarray(X)
        self.X_fit_ = X
        n_samples = X.shape[0]

        # Compute pairwise squared Euclidean distances
        D = squareform(pdist(X, metric="euclidean"))
        D_squared = D**2

        # Double centering
        H = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
        B = -0.5 * H.dot(D_squared).dot(H)

        # Eigen decomposition
        eigvals, eigvecs = np.linalg.eigh(B)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx][: self.n_components]
        eigvecs = eigvecs[:, idx][:, : self.n_components]

        # Store components for out-of-sample projection
        self.L_sqrt_ = np.diag(np.sqrt(eigvals))
        self.eigvecs_ = eigvecs

        # Final embedding
        self.embedding_ = eigvecs.dot(self.L_sqrt_)
        return self

    def fit_transform(self, X):
        """
        Fit the model from data in X and return the embedded coordinates.
        """
        return self.fit(X).embedding_

    def transform(self, X_new):
        """
        The transform method maps new points into the MDS space by computing their distances
        to the training data, centering these distances, and projecting them using the
        training eigenvectors and eigenvalues, approximating their low-dimensional positions.
        """
        if self.embedding_ is None:
            raise RuntimeError("Must call fit or fit_transform before transform.")

        X_new = np.asarray(X_new)
        n = self.X_fit_.shape[0]

        # Compute pairwise distances from new points to training points
        D_new = cdist(X_new, self.X_fit_, metric="euclidean")
        D_sq = D_new**2

        # Use classical MDS trick to approximate the new embedding
        row_mean = D_sq.mean(axis=1)
        col_mean = D_sq.mean(axis=0)
        total_mean = D_sq.mean()

        B_new = -0.5 * (
            D_sq - col_mean[np.newaxis, :] - row_mean[:, np.newaxis] + total_mean
        )

        X_new = B_new.dot(self.eigvecs_).dot(np.linalg.inv(self.L_sqrt_))

        return X_new


# ---------------------------------------------------------------------
# Usage example / demonstration:

if __name__ == "__main__":
    np.random.seed(42)
    X_train = np.random.randn(10, 3)
    X_test = np.random.randn(3, 3)

    mds = CustomMDS(n_components=2)
    X_train_2d = mds.fit_transform(X_train)
    print("Train data embedded shape:", X_train_2d.shape)

    X_test_2d = mds.transform(X_test)
    print("Test data embedded shape: ", X_test_2d.shape)
