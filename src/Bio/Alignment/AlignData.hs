{- | Data structures and helper functions for calculating alignments

   There are two ways to view an alignment: either as a list of edits
   (i.e., insertions, deletions, or substitutions), or as a set of sequences
   with inserted gaps.

   The edit list approach is perhaps more restrictive model but doesn't generalize
   to multiple alignments.

   The gap approach is more general, and probably more commonly used by other software
   (see e.g. the ACE file format).

-}

{-# LANGUAGE ParallelListComp #-}

module Bio.Alignment.AlignData (
    -- * Data types for gap-based alignemnts
       Sequence(..), Gaps, Alignment
    -- * Helper functions
    , extractGaps, insertGaps
    -- * Data types for edit-based alignments
    , Edit(..), EditList, SubstMx, Selector, Chr
    -- * Helper functions
    , columns, eval, isRepl, on
    , showalign, toStrings
    ) where

import qualified Data.ByteString.Lazy as B
import qualified Data.ByteString.Lazy.Char8 as BC
import Bio.Core.Sequence
import Bio.Core.Strand
import Data.List (unfoldr)
import Data.Word
import Data.Char (chr)

-- ----------------------------------------

-- Sequence is a header String, plus SeqData and optional QualData
-- Sequence is an instance of the Bio.Core.Sequence.BioSeq class.

data Sequence = Seq SeqLabel SeqData (Maybe QualData)
                deriving Eq

instance BioSeq Sequence where
  seqid     (Seq lab seq mqual) = SeqLabel $ BC.takeWhile (/= ' ') $ unSL lab
  seqheader (Seq lab seq mqual) = lab
  seqdata   (Seq lab seq mqual) = seq
  seqlength (Seq lab seq mqual) = Offset {unOff = BC.length $ unSD seq}

type Gaps = [Offset]
type Alignment = [(Offset, Strand, Sequence, Gaps)]

-- | Gaps are coded as '*'s, this function removes them, and returns
--   the sequence along with the list of gap positions.
--   note that gaps are positioned relative to the *gapped* sequence 
--   (contrast to stmassembler/Cluster.hs)
extractGaps :: SeqData -> (SeqData, Gaps)
extractGaps str = (SeqData {unSD = (BC.filter (/='*') (unSD str))}, map (\x -> Offset {unOff = x}) $  BC.elemIndices '*' (unSD str))

-- todo: faster to lift concat out of the inner loop?
insertGaps :: Char -> (SeqData, Gaps) -> SeqData
insertGaps c (str', gaps) = SeqData {unSD = go (unSD str') B.empty 0 (map unOff gaps)}
    where go str acc p (next:rest) = let (a,b) = BC.splitAt (next-p) str
                                     in go b (BC.concat [acc,a,BC.pack [c]]) (next+1) rest
          go str acc _ [] = BC.append acc str

-- ----------------------------------------

-- Q&D helper function
showalign a = let (s1,s2) = toStrings a in s1++"\n"++s2

-- | turn an alignment into sequences with '-' representing gaps
-- (for checking, filtering out the '-' characters should return
-- the original sequences, provided '-' isn't part of the sequence
-- alphabet)
toStrings :: EditList -> (String,String)
toStrings [] = ("","")
toStrings (x:xs) = let (a1',a2') = toStrings xs
                       chr' = chr . fromIntegral
                   in case x of Ins c -> ('-':a1', chr' c:a2')
                                Del c -> (chr' c:a1', '-':a2')
                                Repl c1 c2 -> (chr' c1:a1', chr' c2:a2')

-- | The sequence element type, used in alignments.
type Chr = Word8

-- | An Edit is either the insertion, the deletion,
--   or the replacement of a character.
data Edit = Ins Chr | Del Chr | Repl Chr Chr deriving (Show,Eq)

-- | An alignment is a sequence of edits.
type EditList = [Edit]

-- | True if the Edit is a Repl.
isRepl :: Edit -> Bool
isRepl (Repl _ _) = True
isRepl _ = False

-- | A substitution matrix gives scores for replacing a character with another.
--   Typically, it will be symmetric.  It is type-tagged with the alphabet - Nuc or Amino.
type SubstMx t a = (Chr,Chr) -> a

-- | Evaluate an Edit based on SubstMx and gap penalty
eval :: SubstMx t a -> a -> Edit -> a
eval mx g c = case c of Ins _ -> g; Del _ -> g; Repl x y -> mx (x,y)

-- | A Selector consists of a zero element, and a funcition
--   that chooses a possible Edit operation, and generates an updated result.
type Selector a = [(a,Edit)] -> a

-- ------------------------------------------------------------
-- | Calculate a set of columns containing scores
--   This represents the columns of the alignment matrix, but will only require linear space
--   for score calculation.
columns :: Selector a -> a -> Sequence -> Sequence -> [[a]]
columns f z (Seq _ s1 _) (Seq _ s2 _) = columns' f z s1 s2

columns' :: Selector a -> a -> SeqData -> SeqData -> [[a]]
columns' f zero s1 s2 = let
        -- the first column consists of increasing numbers of insertions
        c0 = zero : map (f.return) (zip c0 (map Ins (B.unpack $ unSD s2)))
        -- given the previous column, and the remains of s2, calculate the next column
        mkcol (p0:prev,x) = if B.null x then Nothing
                            else let xi = B.head x
                                     ys = B.unpack $ unSD s2
                                     c  = f [(p0,Del xi)] : [f [del,ins,rep] | del <- zip prev $ repeat (Del xi)
                                                                             | ins <- zip c $ map Ins ys
                                                                             | rep <- zip (p0:prev) $ map (Repl xi) ys]
                                 in Just (c,(c,B.tail x))
    in c0 : unfoldr mkcol (c0, unSD s1)

on c f x y = c (f x) (f y)
