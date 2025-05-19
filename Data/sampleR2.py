# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 12:57:17 2025

@author: Tim
"""

import numpy as np
import matplotlib.pyplot as plt

def generate_grid_point_cloud(step_size):
    
    x_values = np.arange(0, 1+step_size, step_size)
    y_values = np.arange(0, 1+step_size, step_size)
    
    X, Y = np.meshgrid(x_values, y_values)
    
    point_cloud = np.vstack([X.ravel(), Y.ravel()]).T
    
    return point_cloud

def plot_point_cloud(points):
    
    plt.scatter(points[:,0], points[:,1], s = 1)
    plt.title("Dense Point Cloud of the Unit Square")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.xlim(0,1.5)
    plt.ylim(0,1.5)
    plt.show()
    
def save_point_cloud(points, path):
    np.savetxt(path, points, delimiter=',')
    
step_size = 0.02
point_cloud = generate_grid_point_cloud(step_size)
plot_point_cloud(point_cloud)

file_path = '/Users/Tim/Documents/persistent-stratified-main/Data/R2Data.txt'

save_point_cloud(point_cloud, file_path)