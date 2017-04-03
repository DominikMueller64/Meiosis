## Design decisions

Taking decisions on the design and the interface of an
application is often times harder than the actual programing. Here, I outline
the reasons for a few of such decisions. I decided to keep data structures as
simple and "base" (in terms of basic R data structures) as possible. The only
exception is the reference class `Converter`, which I will comment
on later. Of course, I could have used objects for representing the list of
crossover parameters and individuals (in different representations), but there
are probably more disadvantages than advantages here. Usually, users of a
package might dislike if the package forces its own data structurs upon
them and it is often times a hassle, if not a deal-breaker, to convert back and
forth between the different structures. By using simple R lists, conversion to
other formats (for instance, "gametes within chromosomes" or a two-rowed matrix
representing the whole individual) are one-liners. I also (personally) found
nested lists to be a convenient and clean way of representing genomic data. The
reference class `Converter` is an exception, an the reason for its existence is
the unfortunate lack of a proper (hash-) map data structure (like `dict`
[dictionary] in python) that is represented by a SEXP and can be exposed to C++.
For converting from founder alleles (*segmental) to actual alleles
(*genotypic*), we need to retrieve the correct genotypic data of a founder gamete,
which is internally managed by a C++ hashmap. An alternative would be to use a
matrix to represent all genotypic founder data and let the founder alleles be
encoded by the actual indices that will pull the correct row/column from this
matrix. However, I abandoned this possibility because I thought that it would not
fit with other data structures and would impose unnecessary restriction
on the encoding of the founder alleles.
