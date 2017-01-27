##' Factorise columns in a data frame
##'
##' Optionally, can give new labels
##' @param df The data frame (or data table) to manipulate
##' @param labels (optionally), a character vector of new labels for the table elements; the names of the character vector should be (some or all of) the old elements
##' @import data.table
##' @return a data.table with updated columns
factorise_columns <- function(df, labels)
{
    dt <- data.table::copy(data.table::data.table(df))
    column_classes <- vapply(names(dt), function(x) {class(dt[[x]])}, "")
    for (column in colnames(dt)[which(column_classes %in% c("character", "factor"))])
    {
        if (class(dt[[column]]) == "factor")
        {
            old_labels <- levels(dt[, get(column)])
        } else
        {
            old_labels <- unique(dt[, get(column)])
        }
        new_labels <- old_labels
        if (!missing(labels))
        {
            new_labels[which(new_labels %in% names(labels))] <-
                labels[new_labels[which(new_labels %in% names(labels))]]
        }
        dt[, paste(column) := factor(get(column), levels = old_labels,
                                     labels = new_labels)]
    }
    return(dt)
}
