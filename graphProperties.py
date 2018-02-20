import sys

if (len(sys.argv) < 3):
    print('Usage: python graphProperties.py streamFile commFileByCommId');
    sys.exit(0);

# The following two values must be determined from the stream file.
n = 0;
T = 0;

fStream = open(sys.argv[1]);
graphStream = [];

for line in fStream:
	u = map(int, line.split(','));
	graphStream.append(u);
	if (n < u[1]):
		n = u[1];
	if (n < u[2]):
		n = u[2];
	if (T < u[3]):
		T = u[3];
fStream.close();

print('Inferred n = ' + `n` + ', and T = ' + `T`);

g = [];
c = [];

for i in range(n+1):
	g.append([]);
	c.append([]);

for up in graphStream:
	u = up[1];
	v = up[2];
	t = up[3];
	uType = up[0];
	if (uType == 0):	#edge delete
		flag = False;
		for e in g[u]:
			if (e[0] == v) and  (e[2] == -1):
				e[2] = t;
				flag = True;
				break;
		if not flag:
			print('Deleted edge does not exist! '+ str(t));
#			sys.exit(0);
	if (uType == 2):
		for e in g[u]:
			if (e[0] == v) and (e[2] == -1):
				print('Inserted edge already exists!'+`u` + ', '+ `v` + ', ' + `e[1]` + ', ' + `t`);
#				sys.exit(0);
		g[u].append([v,t,-1]);

currComm = -1;
f = open(sys.argv[2]);
for line in f:
	line = line.strip();
	print(line);
	if ':' in line:
		tok = line.split(':');
		currComm = int(tok[0]);
		line = tok[1];
	if ',' in line:
		line = line.lstrip('(');
		line = line.rstrip(')');
		tok = line.split(',');
		u = int(tok[0]);
		s = int(tok[1]);
		e = int(tok[2]);
		c[u].append([currComm,s,e]);

print('Finished reading community file');

t = 0;

while (t<=T):
	print ('For t = ' + `t`);
	gt = [];
	ct = [];
	for i in range(n):
		gt.append([]);
		ct.append([]);
	for i in range(n):
		for e in g[i]:
			if (e[1] <= t and (e[2] >=t or e[2] == -1)):
				src = i;
				dst = e[0];
				gt[src].append(dst);
		for x in c[i]:
			if (x[1] <= t and (x[2] >= t or x[2] == -1)):
				ct[i].append(x[0]);
			if (i == 0) and (t == 0):
				print('(' + ','.join(map(str,x)) + ')');
				print('(' + ','.join(map(str, ct[i])) + ')');
			
	degree = map(lambda x:len(gt[x]),range(n));
	minD = min(degree);
	maxD = max(degree);
	degreeDist = [0 for i in range(maxD-minD+1)];
	for d in degree:
		degreeDist[d-minD] += 1;
	f = open('GraphProps/degreeDist_'+str(t),'w');
	for j in range(minD,maxD+1):
		f.write(`j`+'\t'+`degreeDist[j-minD]`+'\n');
	f.close();
	print('Finished degree distribution');
	commMem = map(lambda x:len(ct[x]),range(n));
	minD = min(commMem);
	maxD = max(commMem);
	commMemDist = [0 for i in range(maxD-minD+1)];
	for d in commMem:
		commMemDist[d-minD] += 1;
	f = open('GraphProps/commMemDist_'+str(t),'w');
	for j in range(minD,maxD+1):
		f.write(`j`+'\t'+`commMemDist[j-minD]`+'\n');
	f.close();
	print('Finished commSize Distribution');
	communities = [];
	numCommunities = 0;
	for i in range(n):
		for j in ct[i]:
			while (numCommunities < (j+1)):
				communities.append(0);
				numCommunities += 1;
			communities[j] += 1;
			#if j in communities.keys():
			#	communities[j] += 1;
			#else:
			#	communities[j] = 1;
	commSizes = map(lambda x: communities[x],range(numCommunities));
	print('numCommunities = ' + `numCommunities`);
	for i in range(numCommunities):
		if (i == 9036):
			print(' '.join(map(str, communities[i])));
		if (commSizes[i]==1):
			print("Community id with size 1 = "+str(i)+", at time = "+str(t));
	minD = min(commSizes);
	maxD = max(commSizes);
	commSizesDist = [0 for i in range(maxD - minD + 1)];
	for d in commSizes:
		commSizesDist[d-minD] += 1;
	f = open('GraphProps/commSizesDist_'+str(t),'w');
	for j in range(minD,maxD+1):
		f.write(`j`+'\t'+`commSizesDist[j-minD]`+'\n');
	f.close();
	t += 250;
