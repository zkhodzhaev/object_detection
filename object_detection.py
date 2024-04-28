# Zulfidin Khodzhaev

import os
import ubermag
import discretisedfield as df
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.transforms as trans
import scipy.stats as stat
import seaborn as sb
import pandas as pd
from sklearn.cluster import DBSCAN

path = "./**"

for path in glob.glob(path, recursive=True):
    if path.endswith("295.out"):
        my_list = path.split('/')[-1]
        path1 = "/home1/08895/bobojon/scratch/correlation"
        
        #check the paths
        print(path)
        print(path1)
        
        files=glob.glob(path + "/m0*.ovf")
        if len(files) < 100:
            continue

        files1 = files.sort()
        files1 = files

        # Function to count skyrmions for a given range of indices
        def count_skyrmions(d2, index_range, center_min, center_max, eps, min_samples):
            data = []
            for i in index_range:
                for j in range(0, 256):
                    if center_min <= d2[i, j] <= center_max:
                        data.append([i, j])

            if not data:  # Check if the data array is empty
                return 0

            data = np.array(data)

            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            dbscan.fit(data)

            return len(set(dbscan.labels_)) - (1 if -1 in dbscan.labels_ else 0)

        # Define the ranges
        leftin = 100
        leftout = 345
        centerin = 345
        centerout = 675
        rightin = 675
        rightout = 920

        center_min = -1
        center_max = -0.2
        eps = 5
        min_samples = 5

        # Initialize lists to store the number of skyrmions for each range
        skyrmions_left = []
        skyrmions_center = []
        skyrmions_right = []

        # Loop through each file
        for file in files1:
            d = df.Field.fromfile(file)
            d1 = d.array
            d2 = d1[0:1024, 0:256, 0, 2]

            n_skyrmions_left = count_skyrmions(d2, range(leftin, leftout), center_min, center_max, eps, min_samples)
            n_skyrmions_center = count_skyrmions(d2, range(centerin, centerout), center_min, center_max, eps, min_samples)
            n_skyrmions_right = count_skyrmions(d2, range(rightin, rightout), center_min, center_max, eps, min_samples)

            skyrmions_left.append(n_skyrmions_left)
            skyrmions_center.append(n_skyrmions_center)
            skyrmions_right.append(n_skyrmions_right)


        # Create a directory to store the data files
        output_directory = path1

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Calculate the total number of skyrmions for each region
        total_skyrmions_left = sum(skyrmions_left)
        total_skyrmions_center = sum(skyrmions_center)
        total_skyrmions_right = sum(skyrmions_right)

        # Specify custom filenames for the data files including the total number of skyrmions and my_list
        left_filename = f"{my_list}_skyrmions_left_total_{total_skyrmions_left}.txt"
        center_filename = f"{my_list}_skyrmions_center_total_{total_skyrmions_center}.txt"
        right_filename = f"{my_list}_skyrmions_right_total_{total_skyrmions_right}.txt"

        # Save the results to files
        with open(os.path.join(output_directory, left_filename), "w") as f:
            for value in skyrmions_left:
                f.write(f"{value}\n")

        with open(os.path.join(output_directory, center_filename), "w") as f:
            for value in skyrmions_center:
                f.write(f"{value}\n")

        with open(os.path.join(output_directory, right_filename), "w") as f:
            for value in skyrmions_right:
                f.write(f"{value}\n")


        left_filepath = os.path.join(output_directory, left_filename)
        center_filepath = os.path.join(output_directory, center_filename)
        right_filepath = os.path.join(output_directory, right_filename)

        with open(left_filepath, "r") as f:
            skyrmions_left = [int(line.strip()) for line in f]

        with open(center_filepath, "r") as f:
            skyrmions_center = [int(line.strip()) for line in f]

        with open(right_filepath, "r") as f:
            skyrmions_right = [int(line.strip()) for line in f]


        # Create the x-axis values (time in ns)
        time = [i for i in range(len(skyrmions_left))]

        # Create horizontally stacked plots
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

        ax1.plot(time, skyrmions_left, label="Left")
        ax1.set_title("Left")
        ax1.set_xlabel("Time (ns)")
        ax1.set_ylabel("Number of Skyrmions")
        ax1.set_xticks(np.arange(0, len(time), len(time) // 10))
        ax1.set_yticks(np.arange(0, max(skyrmions_left + skyrmions_center + skyrmions_right), 2))
        ax1.grid()

        ax2.plot(time, skyrmions_center, label="Center", color="red")
        ax2.set_title("Center")
        ax2.set_xlabel("Time (ns)")
        ax2.set_xticks(np.arange(0, len(time), len(time) // 10))
        ax2.set_yticks(np.arange(0, max(skyrmions_left + skyrmions_center + skyrmions_right), 2))
        ax2.grid()

        ax3.plot(time, skyrmions_right, label="Right", color="green")
        ax3.set_title("Right")
        ax3.set_xlabel("Time (ns)")
        ax3.set_xticks(np.arange(0, len(time), len(time) // 10))
        ax3.set_yticks(np.arange(0, max(skyrmions_left + skyrmions_center + skyrmions_right), 2))
        ax3.grid()
        plt.tight_layout()
        plt.savefig(path1 + "/" + my_list + "_left_center_right_chambers_comparison_dt1ns.png")
        plt.close("all")


