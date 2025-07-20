# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 15:26:03 2025

@author: alexa
"""

import lambertSolver as lt
import dill  # pip install dill

def test_func():
    return lt.lambert_solver([1, 0, 0], [1, 1, 0], 86400, 1.327e11, trajectory='pro')

if __name__ == '__main__':
    print("Trying to pickle lambert_solver...")
    print(dill.pickles(lt.lambert_solver))  # should print True