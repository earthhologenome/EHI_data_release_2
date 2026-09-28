# 2nd EHI Data Release

This is the repository of the 2nd EHI Data Release.

The raw code used for data analysis is in the **Rmd** files stored in the root directory of this repository, while the bookdown-rendered webbook is available at:

[earthhologenome.github.io/EHI_data_release_2](https://earthhologenome.github.io/EHI_data_release_2)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the Rmd documents. Note that the code chunks that require heavy computation have been tuned off using 'eval=FALSE'. To re-render the webbook, you can use the following code:


```r
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```

Metadata is collected from the authenticated `ehi-core` API. Set
`EHI_CORE_EMAIL` and enter the password interactively, or set
`EHI_CORE_EMAIL`/`EHI_CORE_PASSWORD` for a non-interactive render. A current
Firebase ID token can instead be supplied as `EHI_CORE_TOKEN`. See
`01_metadata_access.Rmd` and `R/ehi_core.R`.

When both methods are configured, the collector signs in with the email and
password to avoid using an expired Firebase token.
