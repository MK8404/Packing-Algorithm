#' @title Extract outlier rows based on row half-bandwidth
#'
#' @description
#' Removes rows (and corresponding columns) whose row score is unusually large,
#' where the row score is computed as \code{max(A[[i]] - i)} (a half-bandwidth proxy
#' for row \code{i} in the sparse list representation).
#'
#' @param A Sparse matrix in list form. Supported formats are either
#' \code{list(A = rows)} or \code{rows} directly (as returned by \code{mat2sp}),
#' where each row is an integer vector of nonzero column indices.
#' @param quant Number between 0 and 1. Rows with scores strictly above the
#' \code{1 - quant} quantile are extracted as outliers.
#' @return A list \code{list(B, retained, extracted)} where:
#' \itemize{
#'   \item \code{B}: sparse matrix list after removing extracted rows/columns.
#' Note that the index of \code{B} will be reset to 1, 2, ... for the retained rows.
#'   \item \code{retained}: integer indices of rows kept from the original input.
#'   \item \code{extracted}: integer indices of rows removed as outliers.
#' }
#'
#

extract_out <- function(A, quant = 0.001) {
    if (quant < 0 || quant > 1) {
        stop("quant must be between 0 and 1")
    }
    if (!is.list(A)) {
        stop("A must be a list representing the sparse matrix")
    }

    # Support both formats:
    # 1) list(A = <row-index list>)
    # 2) <row-index list> directly (as returned by mat2sp)
    rows <- if (!is.null(A$A)) A$A else A

    if (!is.list(rows)) {
        stop("A must be a list representing the sparse matrix")
    }

    # Calculate the l1 norm for each row
    l1_norms <- vector("numeric", length(rows))
    sequence <- seq_along(rows)
    for (i in sequence) {
        if (length(rows[[i]]) == 0) {
            l1_norms[i] <- 0
        } else {
            l1_norms[i] <- max(rows[[i]] - i)
        }
    }
    # Determine the threshold for outliers
    threshold <- quantile(l1_norms, probs = 1 - quant)
    # Identify the indices of the rows to be removed
    extracted <- which(l1_norms > threshold)
    retained <- setdiff(sequence, extracted)
    # Create the new matrix B by removing the identified rows and corresponding columns.
    B <- rows[retained]
    for (i in seq_along(B)) {
        B[[i]] <- B[[i]][which(B[[i]] %in% retained)]
    }
    # Reset the indices in B to be 1, 2, ... for the retained rows
    index_map <- setNames(seq_along(retained), retained)
    for (i in seq_along(B)) {
        B[[i]] <- index_map[as.character(B[[i]])]
    }


    return(list(B = B, retained = retained, extracted = extracted))
}
