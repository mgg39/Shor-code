Pkg.add("Gadfly")


using Gadfly


m = 7*3
r = 2
m,r, gcd(m,r)
# >> 21,2,1


mbits = convert(Int,ceil(log2(m)))
# >> 5


auxbits = mbits+5
# >> 10

qbit = zeros(Complex,2^auxbits,2^mbits)


[qbit[nn,powermod(r,nn,m)]= 1/sqrt(2^auxbits) for nn=1:(2^auxbits)]
spy(qbit)


spy(qbit[1:50,:])


col_norm2 = sum(abs2(qbit),1)
plot(x=collect(1:2^mbits), y = col_norm2, Geom.bar)


rmeas = rand()
loc=0
for nn=1:2^mbits
     if rmeas < col_norm2[nn]
          loc = nn; break;
     else
          rmeas -= col_norm2[nn]
     end
end


plot(x=collect(1:100), y = abs(rem_col[1:100]), Geom.bar)


rem_fft=fft(qbit[:,loc])
rem_fft./=sum(abs(rem_col))
plot(x=collect([2:2^auxbits]), y = abs(rem_fft[2:2^auxbits]),Geom.bar)


loc=0
for nn=1:2^auxbits
  if rmeas < abs(rem_col[nn])
      loc = nn; break;
  else
      rmeas -= abs(rem_col[nn])
  end
end

c = loc-1


offset = ceil(m-3*sqrt(m))
diff = floor(m +1 - 2*sqrt(m)) - offset
dist = zeros(convert(Int,diff+1),1)
[dist[nn+1] = abs( (nn+offset)*c/2^auxbits - round((nn+offset)*c/2^auxbits)) for nn=0:diff]


s = indmin(dist) + offset -1
agcd=gcd(convert(Int,s),m)
s,agcd , m/agcd
(9.0,3,7.0)
