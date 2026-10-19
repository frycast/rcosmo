# Suggested DESCRIPTION references

For maintainer review only. No citations have been added to DESCRIPTION.
This file is excluded from the source package by .Rbuildignore.

The reviewer requested author-year citations followed by a DOI, ISBN or URL
in angle brackets, with no space after doi: or https:.

## Recommended: the rcosmo paper

Daniel Fryer, Ming Li and Andriy Olenko (2020).
"rcosmo: R Package for Analysis of Spherical, HEALPix and Cosmological Data".
The R Journal, 12(1), 206-225.

DOI: 10.32614/RJ-2020-012.

Verified sources (publisher article and PDF):
https://journal.r-project.org/articles/RJ-2020-012/
https://journal.r-project.org/articles/RJ-2020-012/RJ-2020-012.pdf

This paper describes the package itself, its HEALPix data structures,
manipulation, visualisation and spatial statistical methods. It is already
cited in README.md. This is the most direct single reference for the
expanded Description.

Suggested sentence:
```text
    The package and its methods are described by Fryer, Li and Olenko (2020)
    <doi:10.32614/RJ-2020-012>.
```

The DOI, authors, year, volume, issue and page range were checked against
the publisher's article metadata; the PDF confirms the title and authors.

## Optional: the HEALPix framework

K. M. Gorski, E. Hivon, A. J. Banday, B. D. Wandelt, F. K. Hansen,
M. Reinecke and M. Bartelmann (2005).
"HEALPix: A Framework for High-Resolution Discretization and Fast Analysis
of Data Distributed on the Sphere".
The Astrophysical Journal, 622(2), 759-771.
DOI: 10.1086/427976.

Verified author-deposited paper record, including the published year,
journal reference and DOI:
https://arxiv.org/abs/astro-ph/0409513

Suggested sentence:
```text
    Uses the HEALPix framework of Gorski et al. (2005)
    <doi:10.1086/427976>.
```

This supports the pixelisation and indexing framework, rather than every
statistical method in rcosmo. The arXiv preprint is dated 2004; use 2005
when citing the published journal article.

## Before submission

Choose whether to add the package paper alone or both references.
Add the chosen sentence(s) to the existing single-paragraph Description,
then replace the pending-references sentence in cran-comments.md with an
accurate statement of the change. Neither proposed citation is applied by
this commit.

CRAN's communication and citation guidance:
https://contributor.r-project.org/cran-cookbook/general_issues.html#communicating-with-cran
https://contributor.r-project.org/cran-cookbook/description_issues.html#references
https://cran.r-project.org/web/packages/policies.html#Re_002dsubmission
