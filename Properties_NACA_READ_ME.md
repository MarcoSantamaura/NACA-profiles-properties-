# READ ME -- NACA profile properties 
Main purpose:  
Read the coordinates of a  NACA airfoils (simmetrical and not) profile from a data file (.dat o .csv) --> which typically contains the combined (x,y) coordinates for the upper and lower edges. The objective is to determine fundamental geometric properties, specially:  
- Area (the profile's surface area)
- Centroid 
- Maximun thickness and its associated x value 
- Maximun camber
- Moment of inertia
- Plot result and draw the profile    

 > Instructions before starting the code:  
    --> dowload from google the dataset of profile (ex.[AirfoilTools](http://airfoiltools.com/airfoil/naca4digit))  
    --> save it in the same folder of code

    
## Functions:
+ <code>selezione_file_finestra()</code>: opens a pop-up window to choose the NACA file --> that gives you <code>percorso</code> the file's path
+ <code>leggi_dati(percorso)</code>: read and setting the data, creating a new array <code>coordinates</code> (x,y).
+ <code>Gauss_Greem_area_MomentiInerzia(coordinates)</code>: provides area, centroid and moments of inerzia. I used Gauss Green's Theorem. 
+ <code>TE_LE_camber(coordinates)</code>: finds leading edge e trailing edge. It calculates camber and thickness, in addition to  createing a graph with the points obtained from linspace for greater accuracy.
+ The last part concerns starting the functions and displaying the results obtained on the screen.


