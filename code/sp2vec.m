function rv=sp2vec(d)
if issparse(d),
	[i j s] = find(d);
	if d(end,end)==0, i = [i; size(d,1)]; j = [j; size(d,2)]; s = [s; 0]; end
	rv = [i j s];
elseif issptensor(d),
	if d(d.size)==0,
		rv = [d.subs d.vals; d.size zeros];
	else
		rv = [d.subs d.vals];
	end
else
	rv=d;
end
end