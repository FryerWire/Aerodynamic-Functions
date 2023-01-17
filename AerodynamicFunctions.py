
import numpy as np

class AerodynamicFunctions():
    """
    Methods
        - stream_function : Generates the stream-function
        - velocity        : 2D source velocity vector that computes the
                            scalar components of u and v
    """

    def __init__(self, strength, x_coords, y_coords, X_meshed, Y_meshed):
        """
        Attributes
            - strength (float) : Strength of the source/sink
            - x_coords (float) : x-coords of the source/sink
            - y_coords (float) : y-coords of the source/sink
            - X_meshed (float) :  2D Numpy array with the x-coords as a meshpoint
            - Y_meshed (float) :  2D Numpy array with the y-coords as a meshpoint
        """

        self.strength = strength
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.X_meshed = X_meshed
        self.Y_meshed = Y_meshed
        
    def stream_function(self, elementary_flow):
        """
        Parameter
            - elementary_flow (string) : source, sink, doublet, or vortex

        Returns
            - Phi (float) : stream function
        """

        elementary_flow = elementary_flow.lower()
        if (elementary_flow == "source") or (elementary_flow == "sink"):
            return (
                (self.strength / (2 * np.pi)) * 
                np.arctan2(
                    (self.Y_meshed - self.y_coords),
                    (self.X_meshed - self.x_coords)
                )
            )   

        elif (elementary_flow == "doublet"):
            return (
                (-self.strength / (2 * np.pi)) * 
                ((self.Y_meshed - self.y_coords) / 
                (((self.X_meshed - self.x_coords) ** 2) + 
                ((self.Y_meshed - self.y_coords) ** 2)))
            )

        else:
            raise Exception("Wrong input, enter one of the following: source, sink, doublet, vortext")

    def velocity(self, elementary_flow):
        """
        Parameter
            - elementary_flow (string) : Source, sink, doublet, vortex

        Returns
            - u (float) : x-component of the doublet velocity vector field
            - y (float) : y-component of the doublet velocity vector field
        """

        elementary_flow = elementary_flow.lower()
        if (elementary_flow == "source") or (elementary_flow == "sink"):
            common_term = ((self.X_meshed - self.x_coords) ** 2) + ((self.Y_meshed - self.y_coords) ** 2)
            u = (self.strength / (2 * np.pi)) * (self.X_meshed - self.x_coords) / common_term
            v = (self.strength / (2 * np.pi)) * (self.Y_meshed - self.y_coords) / common_term

        elif (elementary_flow == "doublet"):
            common_term = ((self.X_meshed - self.x_coords) ** 2) + ((self.Y_meshed - self.y_coords) ** 2)

            u = (
                (-self.strength / (2 * np.pi)) * ((((self.X_meshed - self.x_coords) ** 2) - 
                ((self.Y_meshed - self.y_coords) ** 2)) / (common_term ** 2))
            )

            v = (
                (-self.strength / (2 * np.pi)) * ((2 * (self.X_meshed - self.x_coords) * 
                (self.Y_meshed - self.y_coords)) / (common_term ** 2))
            )

        else:
            raise Exception("Wrong input, enter one of the following: source, sink, doublet, vortext")  

        return u, v      

