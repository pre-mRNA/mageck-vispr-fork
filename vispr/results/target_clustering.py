import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster

from vispr.results.common import templates


class TargetClustering:
    def __init__(self, target_results, topn=1000):
        self.df = pd.DataFrame()
        for condition, selections in target_results.items():
            targets = pd.concat([targets[:]["score"]
                                 for targets in selections.values()])
            self.df[condition] = targets
        self.conditions = self.df.columns
        variance = self.df.var(axis=1)
        variance.sort(ascending=False)
        self.df = self.df.ix[variance[:topn].index]
        self.linkage = linkage(self.df, method="complete")

    def plot_clustering(self, k):
        """Plot k flat clusters."""
        clustering = pd.DataFrame(fcluster(self.linkage, k,
                                           criterion="maxclust"),
                                  columns=["cluster"])
        clustering.index = self.df.index[clustering.index]
        data = pd.merge(self.df, clustering, left_index=True, right_index=True)
        data.sort("cluster", inplace=True)
        data["target"] = data.index

        conditions = []
        for condition in self.df.columns:
            d = data[["cluster", "target", condition]]
            d.columns = ["cluster", "target", "beta"]
            d["condition"] = condition
            conditions.append(d)
        print(conditions)
        data = pd.concat(conditions)
        plt = templates.get_template("plots/target_clustering.json").render(
            data=data.to_json(orient="records"))
        return plt
