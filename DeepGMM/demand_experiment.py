#!/usr/bin/env python
# coding: utf-8

## Adapting DeepGMM for demand experiment of KIV, DIV and DualIV

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import torch
from scenarios.abstract_scenario import AbstractScenario
from methods.demand_model_selection_method import DemandModelSelectionMethod

import json
import pathlib
import time

# Create demand experiment
class Demand(AbstractScenario):
    def __init__(self):
        AbstractScenario.__init__(self)

pathlib.Path("results").mkdir(parents=True, exist_ok=True)

# Run experiment
num_reps = 20
num_test = 2800
num_trains = [25, 500]
rho_vec = [0.1, 0.25, 0.5, 0.75, 0.9]
# Requires minimum 1001 epochs as psi_eval_burn_in = 50 and learning_eval has eval_freq = 20
num_epochs = [3000]
var = [0.1]

results_col = {}
time_col = {}

with torch.cuda.device(0):
    for v in var:
        for e in num_epochs:
            for i, rho in enumerate(rho_vec):
                results = np.zeros((len(num_trains), num_reps))
                time_col[rho] = {}
                for j, num_train in enumerate(num_trains):
                    num_dev = num_train
                    start = time.time()
                    for rep in range(num_reps):
                        print(
                            "rho: {}, num_train: {}, rep: {}".format(
                                rho, num_train, rep
                            )
                        )
                        # Setup
                        scenario = Demand()
                        scenario.setup_c(
                            num_train,
                            num_dev=num_dev,
                            num_test=num_test,
                            seed=rep,
                            rho=rho,
                        )
                        scenario.to_tensor()
                        if torch.cuda.is_available():
                            scenario.to_cuda()
                        train = scenario.get_dataset("train")
                        dev = scenario.get_dataset("dev")
                        test = scenario.get_dataset("test")

                        # Set random seed
                        seed = rep
                        torch.manual_seed(seed)
                        np.random.seed(seed)

                        scenariotrain = scenario.get_dataset

                        method = DemandModelSelectionMethod(
                            num_epochs=e,
                            var=v,
                            f_lr_mult=5.0,
                            enable_cuda=torch.cuda.is_available(),
                        )
                        method.fit(
                            e,
                            train.x, train.z, train.y,
                            dev.x, dev.z, dev.y,
                            g_dev=dev.g,
                            verbose=True,
                        )
                        g_pred_test = method.predict(test.x)
                        mse = float(((g_pred_test - test.g) ** 2).mean())

                        print("---------------")
                        print("Finished {}".format(rho, num_train))
                        print("Test MSE:", np.log10(mse))
                        print("")

                        results[j, rep] = np.log10(mse)

                    time_col[rho][num_train] = time.time() - start

                results_col[rho] = results.tolist()

                fig, ax = plt.subplots()
                ax.boxplot(results.T)
                ax.set_xlabel("Sample Size")
                ax.set_ylabel("log10 MSE")
                ax.set_title("DeepGMM (rho=%.2f)" % rho)
                ax.set_xticklabels(num_trains)
                plt.tight_layout()
                fig.savefig(
                    "results/init_var_{}_epochs_{}_deepgmm_rho_{:.2f}.pdf".format(
                        v, e, rho
                    ),
                    dpi=100,
                )

            with open(
                "results/init_var_{}_epochs_{}_deepgmm_results.json".format(v, e), "w"
            ) as output:
                json.dump(results_col, output)

            with open(
                "results/init_var_{}_epochs_{}_deepgmm_time.json".format(v, e), "w"
            ) as output:
                json.dump(time_col, output)
