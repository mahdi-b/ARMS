from src.ARMS.util.updateGroups import update_groups
a=["1.groups"]
b=["2.groups"]
c=["3.groups"]
update_groups(a,b,".","1_2_rslt")
update_groups(["./1_2_rslt_updated.groups"],c,".","2_3_rslt")

