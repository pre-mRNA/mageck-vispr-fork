import pandas as pd


class TargetClustering:
    def __init__(self, target_results, topn=1000):
        self.df = pd.DataFrame()
        for condition, selections in target_results.items():
            targets = pd.concat([targets[:]["score"] for targets in selections.values()])
            self.df[condition] = targets
        variance = self.df.var(axis=1)
        variance.sort(ascending=False)
        self.df = self.df.ix[variance[:topn].index]
