from classes.ChewbaccaCommand import *
from util.Ungap_Program_Chewbacca import Ungap_Program_Chewbacca


class Ungap_Command(ChewbaccaCommand):

    supported_programs = [Ungap_Program_Chewbacca]
    default_program = Ungap_Program_Chewbacca
    command_name = "Ungap"
