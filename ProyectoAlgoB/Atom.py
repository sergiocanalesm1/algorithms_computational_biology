#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:29:35 2020

@author: David Ricardo Figueora Blanco 
@email: dr.figueroa10@uniandes.edu.co 

   FUNCTION:  a short paragraph stating the purpose of the
              program.

   INPUT:     location of input, i.e.  standard input, a file on
              disk

   OUTPUT:    location and type of output, i.e.  a report
              containing a detail record for each city processed
              containing city id, Celsius temperature, Fahrenheit
              temperature and wind chill temperature.

   NOTES:     any relevant information that would be of
              additional help to someone looking at the program

"""
## modules
import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys

class Atom:

    def __init__(self,atom_name,atom_type,atom_id,res_name,res_num,coords):
        """
        Parameters
        ----------
        name : string
            DESCRIPTION.
        res_name : string
            DESCRIPTION.
        res_num : int
            DESCRIPTION.
        coords : list
            [x,y,z] coordiantes of atom
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.atom_name = atom_name
        self.atom_type = atom_type
        self.res_name = res_name
        self.res_num = res_num
        self.coords = coords
        self.atom_id = atom_id


    def get_Name(self):
        """
        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.atom_name
    
    def get_atom_type(self):
        return self.atom_type
    
    def get_resName(self):
        """


        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.res_name
    def get_resNum(self):
        """


        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.res_num
    def get_Coords(self):
        """


        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.coords
    def get_Atomid(self):
        return self.atom_id

