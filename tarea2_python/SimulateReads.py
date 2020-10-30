#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:24:39 2020

@author: David Ricardo Figueora Blanco 
@email: dr.figueroa10@uniandes.edu.co 

   FUNCTION:  Simulate reads from fasta file
   INPUT:     location of input, i.e.  standard input, a file on
              disk

   OUTPUT:    location and type of output, i.e.  a report
              containing a detail record for each city processed
              containing city id, Celsius temperature, Fahrenheit
              temperature and wind chill temperature.

   NOTES:     any relevant information that would be of
              additional help to someone looking at the program

"""
import random


    
def SimulateReads(seq="ATGGTACTGACGATTTATCCT",readlen=5,reads=9):
    simulated_list=[]
    SeqLen=len(seq)

    
    for i in range(reads):
        rand_num=random.randint(0,SeqLen-readlen)
        #print(seq[rand_num:rand_num+readlen])
        simulated_list.append(seq[rand_num:rand_num+readlen])
    return simulated_list
