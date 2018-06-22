function sc(x)
  if(abs(x)>1.0)
	a=trunc(log10(abs(x)))
	b=x/(10.0^a)
  else
	a=trunc(log10(abs(x)))-1
	b=x/(10.0^a)
  end
  return a,b
end
