
a = ""


還有分成正反股

contextdict={}

if a[1]=="G":
	contextdict[a]="CG"

if a[1]!="G":
	if a[2]=="G":
		contextdict[a]="CHG"
	if a[2]!="G":
		contextdict[a]="CHH"
