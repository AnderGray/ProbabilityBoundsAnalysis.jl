# Arithmetic 

<!---

Supported dependent arithmetic between uncertain numbers:

|                           |     independent    | dependency known   | dependency unknown | perfect/opposite     | partial information  |
|:-------------------------:|:------------------:|:------------------:|:------------------:|:--------------------:|:--------------------:|
| intervals                 | not known to exist |         solutions exist        |         yes        |  solutions exist |  solutions exist |
| probability distributions |         yes        |         yes        |         yes        |          yes         |          yes        |
| probability boxes         |         yes        |         yes        |         yes        |          yes         |          yes        |


||independent|dependency known|dependency unknown|perfect/opposite|partial information|
|---------------------------|:------------------:|--------------------|--------------------|----------------------|----------------------|
|intervals|not known to exist|solutions exist|yes|solutions exist|solutions exist|
|probability distributions|yes|yes|yes|yes|yes|
|probability boxes|yes|yes|yes|yes|yes|

--->

Most of the fundamental binary operations can be performed between uncertain numbers of all types:

```julia
julia> a = normal(-1,1); 
julia> b = interval(1,2);
julia> a + b
 Pbox: 	  ~  ( range=[-3.09023,4.090232], mean=[0.0,1.0], var=[1.0,1.25])

julia> a - b
 Pbox: 	  ~  ( range=[-6.090232,1.0902323], mean=[-3.0,-2.0], var=[1.0,1.25])

julia> a * b
 Pbox: 	  ~  ( range=[-4.090232,4.180464], mean=[-1.015451,-1.9690], var=[0.99763,3.99053])

julia> a / b
 Pbox: 	  ~  ( range=[-2.045116,2.0902323], mean=[-0.50772,-0.984548], var=[0.249408,0.99763])
```

All of the above operations assume independence. For unknown dependence:
```julia
julia> convFrechet(a, b, op = +)
 Pbox: 	  ~  ( range=[-3.09023,4.090232], mean=[0.0,1.0], var=[0.383917,2.1086384])

julia> convFrechet(a, b, op = -)
 Pbox: 	  ~  ( range=[-6.09023,1.090232], mean=[-3.0,-2.0], var=[0.383917,2.108638])
```

The resulting p-boxes are much wider than the independence case.

Perfect and opposite convolutions can also be performed:
```julia
julia> a = normal(0,1);
julia> b = normal(1,1);
julia> convPerfect(a, b, op = +)
 Pbox: 	  ~  ( range=[-5.18046,7.18046], mean=[0.96909,1.030903], var=[3.80050,4.18248])

julia> convOpposite(a, b, op = +)
 Pbox: 	  ~  ( range=[0.48559,1.51440], mean=[0.96909,1.03090], var=[0.0,0.00840])
```

Binary operations with a specified correlation coefficient may also be performed:

```julia
julia> a = normal(0,1);
julia> b = normal(1,1);
julia> conv(a,b, op = +, corr = 0.5)
 Pbox: 	  ~  ( range=[-5.18046,7.18046], mean=1.0, var=[2.57835,3.96457])
```
This assumes that a and b follow a Gaussian Copula. You may however perform the operation with any copula by using the function
```julia
julia> convCorr(a, b, C = C, op = +)
```
where C is a copula (see section on dependence modelling).

Note that:

```julia
conv(a, b, op = op, corr = 0)               == convIndep(a, b, op = op)
conv(a, b, op = op, corr = 1)               == convPerfect(a, b, op = op)
conv(a, b, op = op, corr = -1)              == convOpposite(a, b, op = op)
conv(a, b, op = op, corr = interval(-1,1))  == convFrechet(a, b, op = op)
```
