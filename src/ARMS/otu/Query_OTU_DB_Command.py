from classes.ChewbaccaCommand import *
from otu.Query_OTU_DB_Program_Vsearch import Query_OTU_DB_Program_Vsearch


class Query_OTU_DB_Command(ChewbaccaCommand):

    supported_programs = [Query_OTU_DB_Program_Vsearch]
    default_program = Query_OTU_DB_Program_Vsearch
    command_name = "Query OTU Identity"
