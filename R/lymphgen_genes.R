#' @title Genes whose mutations a LymphGen model can use
#'
#' @description The gene symbols, as they appear in a MAF's Hugo_Symbol
#' column, of every gene a LymphGen model has a feature for, so a MAF can be
#' restricted to them before [lymphgen_input], e.g. while it is read with
#' `GAMBLR.results::get_all_coding_ssm(..., these_genes = lymphgen_genes())`.
#' Mutations of other genes can't affect the model's calls.
#'
#' @details The genes are those of every feature in the model's training
#' feature annotation (`Full.Uber.2020.txt` in the LymphGen repository, or in
#' `models/<model>/` for the GAMBL models), copy number features included,
#' plus the genes of the LymphGen features of the GAMBLR.utils classifier
#' registry, which [lymphgen_input] also uses. Pseudo-genes are given as the
#' gene their mutations come from: MYD88 for MYD88 L265P and non-L265P,
#' NFKBIZ for the NFKBIZ 3' UTR, FOXO1 for FOXO1 N-terminal, and the genes in
#' `split_genes.txt` for a model's split pseudo-genes. A pseudo-gene with no
#' known parent gene is an error rather than being left out, since leaving it
#' out would silently drop its mutations.
#'
#' Each gene is returned under every symbol [lymphgen_input] maps to its
#' Entrez ID (the registry's symbols and those of `GAMBLR.utils::gene_tss`),
#' so a MAF filtered with these symbols keeps every row [lymphgen_input]
#' would use. The gene list LymphGen refits its model to doesn't come from
#' the MAF, so filtering the MAF doesn't change it.
#'
#' @param model "original" (default), "gambl_v2", "gambl_v3", "gambl_v4" or
#' "GAMBL-U".
#' @param lymphgen_dir The LymphGen repository. Default: the
#' `GAMBLR.predict.lymphgen_dir` option, else `LYMPHGEN_DIR` (see
#' [run_lymphgen]).
#'
#' @return A sorted character vector of gene symbols.
#' @export
#'
#' @examples
#' \dontrun{
#' options(GAMBLR.predict.lymphgen_dir = "~/git/Lymphgen_2")
#' meta = get_gambl_metadata() %>%
#'   dplyr::filter(pathology == "DLBCL", seq_type %in% c("genome", "capture"))
#' # read only the mutations GAMBL-U can use
#' maf = get_all_coding_ssm(meta, include_silent = TRUE,
#'                          these_genes = lymphgen_genes("GAMBL-U"))
#' }
lymphgen_genes = function(model = "original", lymphgen_dir = NULL){
  model = match.arg(model, c("original", "gambl_v2", "gambl_v3", "gambl_v4", "GAMBL-U"))
  path = .lymphgen_uber_path(model, lymphgen_dir)
  catalogue = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  ids = unique(as.numeric(catalogue$ENTREZ.ID[catalogue$ENTREZ.ID > 0]))
  # pseudo-genes (Entrez IDs from 1e9) -> the gene whose mutations they are
  parent = c(`1000000004` = 4615,                              # MYD88 L265P
             `1000000005` = 4615,                              # MYD88 non-L265P (.gambl_v2)
             `1000000006` = 64332,                             # NFKBIZ 3' UTR (.gambl_v2)
             `1000000007` = 2308)                              # FOXO1 N-terminal (.gambl_v3)
  split_file = file.path(dirname(path), "split_genes.txt")
  if(model != "original" && file.exists(split_file)){
    split = utils::read.delim(split_file, stringsAsFactors = FALSE)
    parent = c(parent, stats::setNames(as.numeric(split$entrez), as.character(split$pseudo)))
  }
  pseudo = ids[ids >= 1e9]
  unknown = setdiff(format(pseudo, scientific = FALSE, trim = TRUE), names(parent))
  if(length(unknown) > 0){
    stop("lymphgen_genes(): pseudo-gene(s) ", paste(unknown, collapse = ", "), " of model ", model,
         " have no known parent gene; add them to lymphgen_genes()")
  }
  entrez = unique(c(ids[ids < 1e9], parent[format(pseudo, scientific = FALSE, trim = TRUE)]))
  known = .lymphgen_symbols()
  registry = GAMBLR.utils::get_classifier_features("LymphGen")
  sort(unique(c(known$gene[known$entrez %in% entrez], stats::na.omit(registry$gene))))
}

# The symbol -> Entrez ID table lymphgen_input() maps MAF genes with: the
# registry's LymphGen genes first, then GAMBLR.utils::gene_tss (grch37,
# primary), one Entrez ID per symbol. Shared with lymphgen_genes() so the two
# can't disagree.
.lymphgen_symbols = function(){
  registry = GAMBLR.utils::get_classifier_features("LymphGen")
  tss = GAMBLR.utils::gene_tss
  tss = tss[tss$genome_build == "grch37" & tss$primary, , drop = FALSE]
  known = unique(rbind(registry[!is.na(registry$gene) & !is.na(registry$entrez), c("gene", "entrez")],
                       tss[, c("gene", "entrez")]))
  known[!duplicated(known$gene), , drop = FALSE]
}
