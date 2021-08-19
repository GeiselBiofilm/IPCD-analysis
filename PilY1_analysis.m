% load BLAST results
data_o1=blastreadlocal('pilY1_O1.out',0);
data_14=blastreadlocal('pilY1_14.out',0);
data_83=blastreadlocal('pilY1_83.out',0);

% load amino acid sequences
a14=fastaread('pilY1_14_aa.txt');
a14.Sequence=a14.Sequence(a14.Sequence~=' ');
ao1=fastaread('pilY1_O1_aa.txt');
ao1.Sequence=ao1.Sequence(ao1.Sequence~=' ');

% find complete amino acid sequences in PA14 clade
for i=1:131
    s14(i).Name=data_14.Hits(i).Name;
    ss=data_14.Hits(i).HSPs(1).Alignment(3,:);
    ss=ss(find(ss~='-'));
    for j=1:3
        s(j)=swalign(nt2aa(ss(j:end),'ACGTOnly','false'),a14.Sequence);
    end
    [bla,j]=max(s);
    s14(i).Sequence=nt2aa(ss(j:end),'ACGTOnly','false');
    l14(i)=length(ss);
end
s14=s14(find(l14>3000));

% find complete amino acid sequences in PAO1 clade
% we use results from both PAO1 and IPCD83 BLAST searches
for i=1:726
    so1(i).Name=data_o1.Hits(i).Name;
    ss=data_o1.Hits(i).HSPs(1).Alignment(3,:);
    ss=ss(find(ss~='-'));
    for j=1:3
        s(j)=swalign(nt2aa(ss(j:end),'ACGTOnly','false'),ao1.Sequence);
    end
    [bla,j]=max(s);
    so1(i).Sequence=nt2aa(ss(j:end),'ACGTOnly','false');
    lo1(i)=length(ss);
end
so1=so1(find(lo1>3000));

for i=1:726
    s83(i).Name=data_83.Hits(i).Name;
    ss=data_83.Hits(i).HSPs(1).Alignment(3,:);
    ss=ss(find(ss~='-'));
    for j=1:3
        s(j)=swalign(nt2aa(ss(j:end),'ACGTOnly','false'),ao1.Sequence);
    end
    [bla,j]=max(s);
    s83(i).Sequence=nt2aa(ss(j:end),'ACGTOnly','false');
    l83(i)=length(ss);
end
s83=s83(find(l83>3000));

% calculate Shannon diversity in PA14 clade
pily14.align=multialign([s14]);
pily14.profile=seqprofile(pily14.align,'Alphabet','AA','Gaps','all');
pily14.shannon=mutated(pily14.profile);

% calculate Shannon diversity in PAO1 clade
stot=[so1,s83];
pilyo1.align=multialign([stot]);
pilyo1.profile=seqprofile(pilyo1.align,'Alphabet','AA','Gaps','all');
pilyo1.shannon=mutated(pilyo1.profile);

% align all PilY1 sequences
sall=[s14,so1,s83];
align_all=multialign([sall]);

% generate phylogenetic tree of PilY1 sequences
D = seqpdist(align_all);
for i=1:852
    names{i}=num2str(i);
end
names{1}='PA14';
names{130}='PAO1';
names{168}='IPCD83';
tree=seqlinkage(D,'average',names);
phytreewrite('tree.tree', tree)


% function to calculate Shannon diversity
function m=mutated(M)
    s=@(v) sum(-v(v>0).*log(v(v>0)));
    for i=1:size(M,2)
        m(i)=s(M(:,i));
    end
end
