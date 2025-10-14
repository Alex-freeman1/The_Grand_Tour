# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 13:19:45 2025

@author: alexa
"""

#Date calculator

from datetime import datetime, timedelta

# # Ask for input date in YYYY-MM-DD format
# input_date_str_dep = "1976-02-01"
# input_date_str_arr = "1977-08-20"


# # Convert string to datetime object
# input_date_dep = datetime.strptime(input_date_str_dep, "%Y-%m-%d")
# input_date_arr = datetime.strptime(input_date_str_arr, "%Y-%m-%d")


# # Define time delta (200 days)
# del_t_dep = 2025
# del_t_arr = 2640
# delta_dep = timedelta(days=del_t_dep)
# delta_arr = timedelta(days=del_t_arr)

# # Calculate new date
# new_date_dep = input_date_dep + delta_dep
# new_date_arr = input_date_arr + delta_arr

# # Print result
# print(f"The date {del_t_dep} days after", input_date_dep.strftime("%Y-%m-%d"), "is", new_date_dep.strftime("%Y-%m-%d"))
# print(f"The date {del_t_arr} days after", input_date_arr.strftime("%Y-%m-%d"), "is", new_date_arr.strftime("%Y-%m-%d"))



date1_str = "1983-01-01"
date2_str = "1989-08-25"


#Convert strings to datetime objects
date1 = datetime.strptime(date1_str, "%Y-%m-%d")
date2 = datetime.strptime(date2_str, "%Y-%m-%d")

# Compute the difference in days
difference = (date1 - date2).days

print(f"The number of days between is {abs(difference)} days.")


