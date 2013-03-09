import re,sys,string,collections

f = open('parasites.txt', 'r')
all = f.read()
f.close()

tmp = re.split("\r\r\r ", all)
res = re.findall('(Parasite: )([\w,\s]+)', tmp[0])


def readFile(filename):
	f = open(filename, 'r')
	alldat = f.read()
	dat = re.split("\r\r\r ", alldat)
	f.close()
	return(dat)

def parasiteLine(record):
	info = re.findall('Parasite: (\w+ \w+[ \w+]?)\W+(\w+, \d+)', record)
	return(info)

def fromLine(record):
	info = re.findall('From: (\w+ \w+[ \w+]?)\W+([\w\s]+)', record)
	return(info)

def localityLine(record):
	info = re.findall('Locality: ([\w\s]+)', record)
	return(info)

def commentLine(record):
	info = re.findall('Comment: ([\w\s]+)', record)
	return(info)

def refLine(record):
	info = re.findall('Ref.No.\W+(\d+){1}\W+([\w\s]+)', record)
	return(info)

refLine(tmp[0])


def parseFile(dat):
	parasites  = [0]*len(dat)
	prefs      = [0]*len(dat)
	hosts      = [0]*len(dat)
	sources    = [0]*len(dat)
	localities = [0]*len(dat)
	comments   = [0]*len(dat)
	refno      = [0]*len(dat)
	refinfo    = [0]*len(dat)
	for i in range(0,len(dat)):
		temp = parasiteLine(dat[i])
		parasites[i] = temp[0][0]
		prefs[i] = temp[0][1]
		temp = fromLine(dat[i])
		hosts[i] = temp[0][1]
		sources[i] = temp[0][1]
		temp = localityLine(dat[i])
		localities[i] = temp
		temp = commentLine(dat[i])
		comments[i] = temp
		temp = refLine(dat[i])
		refno[i] = temp[0][0]
		refinfo[i] = temp[0][1]
	output = collections.namedtuple('group', ['parasites', 'prefs', 'hosts', 'sources', 'localities', 'comments', 'refno', 'refinfo'])
	a = output(parasites, prefs, hosts, sources, localities, comments, refno, refinfo)
	return(a)


def main():
	filename = sys.argv[1]
	print filename
	dat = readFile(filename)
	p = parseFile(dat)
	print p

if __name__ == "__main__":
	main()
