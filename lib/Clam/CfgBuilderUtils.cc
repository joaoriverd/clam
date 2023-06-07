#include "CfgBuilderUtils.hh"

#include "llvm/ADT/APInt.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/PostOrderIterator.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/IR/CallSite.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/IR/DebugInfo.h"
#include "llvm/IR/DebugLoc.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/Value.h"

#include <cstdint>
#include <algorithm>
#include <set>

namespace clam {

using namespace ikos;
using namespace llvm;

// Any integer that cannot be represented by 64 bits is considered a bignum.
bool isSignedBigNum(const APInt &v) {
  unsigned b = v.getBitWidth();
  if (b <= 64) {
    return false;
  } else {
    // if bitwidth > 64 then we check the actual value
    APInt max(b, APInt::getSignedMaxValue(64).getSExtValue(), true);
    APInt min(b, APInt::getSignedMinValue(64).getSExtValue(), true);
    return (v.sgt(max) || v.slt(min));
  }
}

bool isBool(const Type *t) { return (t->isIntegerTy(1)); }

bool isBool(const Value &v) { return isBool(v.getType()); }

bool isInteger(const Type *t) { return (t->isIntegerTy() && !isBool(t)); }

bool isInteger(const Value &v) { return isInteger(v.getType()); }

bool isDouble(const Type *t) { return t->isDoubleTy(); }

bool isDouble(const Value &v) { return isDouble(v.getType()); }

bool isFloat(const Type *t) { return t->isFloatTy(); }

bool isFloat(const Value &v) { return isFloat(v.getType()); }

bool isFloatingPoint(const Type *t) { return t->isFloatingPointTy(); }

bool isFloatingPoint(const Value &v) { return isFloatingPoint(v.getType()); }

bool isReference(const Type *t, const CrabBuilderParams &params) {
  return (t->isPointerTy() && params.trackMemory());
}

bool isReference(const Value &v, const CrabBuilderParams &params) {
  return isReference(v.getType(), params);
}

z_number toZNumber(const APInt &v, const CrabBuilderParams &params,
		   bool interpretAsSigned, bool &isTooBig) {                   
  isTooBig = false;
  if (!params.enable_bignums) {
    isTooBig = isSignedBigNum(v);
  }
#if 0
  // Convert to strings is not ideal but it shouldn't be a big
  // bottleneck.
  std::string val = v.toString(10, interpretAsSigned);
  return z_number(val);
#else
  // Based on:
  // https://llvm.org/svn/llvm-project/polly/trunk/lib/Support/GICHelper.cpp  
  if (!interpretAsSigned) {
    return z_number::from_raw_data(v.getRawData(), v.getNumWords());
  } else {    
    APInt abs;
    abs = v.isNegative() ? v.abs() : v;
    const uint64_t *rawdata = abs.getRawData();
    unsigned numWords = abs.getNumWords();
    z_number res(z_number::from_raw_data(rawdata, numWords));
    return v.isNegative() ? -res : res;
  }
#endif
}

z_number toFPNumber(const APFloat &v, const CrabBuilderParams &params) {
    return z_number(v.convertToDouble(), FP_Number);
}

z_number getIntConstant(const ConstantInt *CI, const CrabBuilderParams &params,
			bool interpretAsSigned, bool &isTooBig) {
                        
  isTooBig = false;
  if (CI->getType()->isIntegerTy(1)) {
    return z_number((int64_t)CI->getZExtValue());
  } else {
    return toZNumber(CI->getValue(), params, interpretAsSigned, isTooBig);
  }
}

z_number getFPConstant(const ConstantFP *CI, const CrabBuilderParams &params) {
    return toFPNumber(CI->getValueAPF(), params);
}

bool isTrackedType(const Type &ty, const CrabBuilderParams &params) {
  return isReference(&ty, params) || ty.isIntegerTy() || ty.isDoubleTy() || ty.isFloatTy();
}

bool isTracked(const Value &v, const CrabBuilderParams &params) {
  // -- ignore any shadow variable created by seahorn
  // if (v.getName().startswith("shadow.mem"))
  // return false;

  return isTrackedType(*v.getType(), params);
}

bool ShouldCallSiteReturn(CallInst &I, const CrabBuilderParams &params) {
  CallSite CS(&I);
  if (Function *Callee =
          dyn_cast<Function>(CS.getCalledValue()->stripPointerCasts())) {
    Type *RT = Callee->getReturnType();
    return (!(RT->isVoidTy()) && isTrackedType(*RT, params));
  }
  return false;
}

bool DoesCallSiteReturn(CallInst &I, const CrabBuilderParams &params) {
  return (!I.getType()->isVoidTy() && isTracked(I, params));
}

bool hasDebugLoc(const Instruction *inst) {
  if (!inst)
    return false;
  const DebugLoc &dloc = inst->getDebugLoc();
  return dloc;
}

crab::cfg::debug_info getDebugLoc(const Instruction *I) {
  return getDebugLoc(I, 0);
}

crab::cfg::debug_info getDebugLoc(const Instruction *I, uint32_t Id) {
  if (hasDebugLoc(I)) {
    const DebugLoc &dloc = I->getDebugLoc();
    unsigned Line = dloc.getLine();
    unsigned Col = dloc.getCol();
    std::string File = (*dloc).getFilename();
    return crab::cfg::debug_info(File, Line, Col, Id == 0 ? -1 : (int64_t) Id);
  } else {
    return crab::cfg::debug_info(Id);
  } 
}

uint64_t storageSize(const Type *t, const DataLayout &dl) {
  return dl.getTypeStoreSize(const_cast<Type *>(t));
}

void normalizeCmpInst(CmpInst &I) {
  switch (I.getPredicate()) {
  case ICmpInst::ICMP_UGT:
  case ICmpInst::ICMP_SGT:
  case FCmpInst::FCMP_OGT:
  case FCmpInst::FCMP_UGT:
    I.swapOperands();
    break;
  case ICmpInst::ICMP_UGE:
  case ICmpInst::ICMP_SGE:
  case FCmpInst::FCMP_OGE:
  case FCmpInst::FCMP_UGE:
    I.swapOperands();
    break;
  default:;
  }
}

bool isIntToBool(const CastInst &I) {
  return (isa<TruncInst>(I) && I.getDestTy()->isIntegerTy(1));
}

bool isBoolToInt(const CastInst &I) {
  return ((isa<ZExtInst>(I) || isa<SExtInst>(I)) &&
          I.getSrcTy()->isIntegerTy(1));
}

bool isBoolArray(const Type &T) {
  return (T.isArrayTy() && T.getArrayElementType()->isIntegerTy(1));
}

bool isIntArray(const Type &T) {
  return (T.isArrayTy() && T.getArrayElementType()->isIntegerTy() &&
          !(T.getArrayElementType()->isIntegerTy(1)));
}

// bool isPointerArray(const Type &T) {
//   return (T.isArrayTy() && T.getArrayElementType()->isPointerTy());
// }

bool isAssertFn(const Function &F) {
  return (F.getName().equals("verifier.assert") ||
          F.getName().equals("crab.assert") ||
	  F.getName().equals("__VERIFIER_assert") || 
          F.getName().equals("__CRAB_assert"));
}

bool isSeaHornFail(const Function &F) {
  return (F.getName().equals("seahorn.fail"));
}

bool isErrorFn(const Function &F) {
  return (F.getName().equals("seahorn.error") ||
          F.getName().equals("verifier.error") ||
          F.getName().equals("__VERIFIER_error") ||
          F.getName().equals("__SEAHORN_error"));
}

bool isAssumeFn(const Function &F) {
  return (F.getName().equals("verifier.assume") ||
          F.getName().equals("__VERIFIER_assume") ||
	  F.getName().equals("__SEA_assume") ||
          F.getName().equals("__CRAB_assume") ||
	  F.getName().equals("llvm.assume"));
}

bool isNotAssumeFn(const Function &F) {
  return (F.getName().equals("verifier.assume.not") ||
          F.getName().equals("__VERIFIER_assume_not") ||
          F.getName().equals("__CRAB_assume_not"));
}

bool isVerifierCall(const Function &F) {
  return (F.isDeclaration() && 
	  (isAssertFn(F) || isErrorFn(F) || isAssumeFn(F) || isNotAssumeFn(F) ||
	   isSeaHornFail(F)));
}

static bool isSeaHornIntrinsic(const Function &F) {
  return (F.isDeclaration() &&
	  (F.getName() == "sea.is_dereferenceable" ||
	   F.getName() == "sea_is_dereferenceable"));
}

bool isCrabIntrinsic(const Function &F) {
  return (F.isDeclaration() &&
	  (F.getName().startswith("__CRAB_intrinsic_") ||
	   isSeaHornIntrinsic(F)));
}

static const std::set<std::string> math_functions {
  "fabs", "fmax", "fmin", "exp", "exp2", "expm1", "log", "log2", "log10", "log1p",
  "sqrt", "cbrt", "pow", "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
  "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "llvm.fabs.f64", // Special names for some math functions
  "__CRAB_get_range", // Special function for debugging. To remove later.
  "__CRAB_get_planar_poly", // Special function for debugging. To remove later.
  "__CRAB_disable_accurate_math", // Special function for debugging. To remove later.
  "__CRAB_accuracy_level_1_math", // Special function for debugging. To remove later.
  "__CRAB_disable_sound_rounding", // Special function for debugging. To remove later.
  "__CRAB_enable_sound_rounding", // Special function for debugging. To remove later.
  "__CRAB_get_max_tvpi_error" // Special function for debugging. To remove later.
};

bool isMathFunction(const Function &F) {
  std::string name =  F.getName();
  bool isMathFun = (math_functions.find(name) != math_functions.end());
  return (F.isDeclaration() && isMathFun);
}

std::string getCrabIntrinsicName(const Function &F) {
  assert(isCrabIntrinsic(F));
  
  if (isSeaHornIntrinsic(F)) {
    // if the bitcode has been processed already by SeaHorn
    StringRef res = F.getName().split("sea.").second;
    if (res == "") {
      // if the user writes the intrinsic on the C source
      res = F.getName().split("sea_").second;
    }
    return res.str();
  } else {
    StringRef res = F.getName().split("__CRAB_intrinsic_").second;
    return res.str();
  }
}

bool isZeroInitializer(const Function &F) {
  return F.getName().startswith("verifier.zero_initializer");
}

bool isZeroInitializer(const CallInst &CI) {
  ImmutableCallSite CS(&CI);
  const Value *calleeV = CS.getCalledValue();
  if (const Function *callee =
          dyn_cast<Function>(calleeV->stripPointerCasts())) {
    return isZeroInitializer(*callee);
  }
  return false;
}

bool isIntInitializer(const Function &F) {
  return F.getName().startswith("verifier.int_initializer");
}

bool isIntInitializer(const CallInst &CI) {
  ImmutableCallSite CS(&CI);
  const Value *calleeV = CS.getCalledValue();
  if (const Function *callee =
          dyn_cast<Function>(calleeV->stripPointerCasts())) {
    return isIntInitializer(*callee);
  }
  return false;
}

std::string getAssertKindFromMetadata(MDNode *MDN) {
  // assume MDN is the metadata associate to getMetadata("clam-assertion")
  if (MDN) {
    if (MDTuple *t = dyn_cast<MDTuple>(MDN->getOperand(0))) {
      std::string checkName = cast<MDString>(t->getOperand(0))->getString();
      return checkName;
    } 
  }
  return "";
}

// Return true if any use is a verifier call  
bool AnyUseIsVerifierCall(Value &V) {
  for (auto &U : V.uses()) {
    Value *User = U.getUser();
    if (CallInst *CI = dyn_cast<CallInst>(User)) {
      CallSite CS(CI);
      const Function *callee = dyn_cast<Function>(CS.getCalledValue()->stripPointerCasts());
      if (callee && isVerifierCall(*callee)) {
	return true;
      }
    }
  }
  return false;
}
  
// Return true if all uses are BranchInst's
bool AllUsesAreBrInst(Value &V) {
  // XXX: do not strip pointers here
  for (auto &U : V.uses())
    if (!isa<BranchInst>(U.getUser()))
      return false;
  return true;
}

// Return true if all uses are either BranchInst's or Select's if
// satisfies selectFilter.
bool AllUsesAreBrOrIntSelectCondInst(Value &V,
                                     const CrabBuilderParams &params,
				     std::function<bool(SelectInst*)> selectFilter) {
  // XXX: do not strip pointers here
  for (auto &U : V.uses()) {
    if ((!isa<BranchInst>(U.getUser())) && (!isa<SelectInst>(U.getUser()))) {
      return false;
    }
    if (SelectInst *SI = dyn_cast<SelectInst>(U.getUser())) {
      if (isBool(*SI) || SI->getCondition() != &V || isReference(*SI, params)) {
        // if the operands are bool or V is not the condition
        return false;
      }
      if (!selectFilter(SI)) {
	return false;
      }
    }
  }
  return true;
}

// Return true if all uses are the callee at callsites
bool AllUsesAreIndirectCalls(Value &V) {
  // XXX: do not strip pointers here
  for (auto &U : V.uses()) {
    if (CallInst *CI = dyn_cast<CallInst>(U.getUser())) {
      CallSite CS(CI);
      const Value *callee = CS.getCalledValue();
      if (callee == &V)
        continue;
    }
    return false;
  }
  return true;
}
  
// Return true if all uses are verifier calls (assume/assert)
bool AllUsesAreVerifierCalls(Value &V, bool goThroughIntegerCasts,
                             bool nonBoolCond,
                             SmallVector<CallInst *, 4> &verifierCalls,
			     bool onlyAssume) {

  auto isVeriCall = [&onlyAssume](const Function &F) {
    if (!onlyAssume) {
      return isVerifierCall(F);
    } else {
      return (isAssumeFn(F) || isNotAssumeFn(F));
    }
  };
    
  for (auto &U : V.uses()) {
    Value *User = U.getUser();
    if (goThroughIntegerCasts) {
      if (isa<ZExtInst>(User) || isa<SExtInst>(User)) {
        return AllUsesAreVerifierCalls(*User, goThroughIntegerCasts,
                                       nonBoolCond, verifierCalls, onlyAssume);
      }
    }

    if (CallInst *CI = dyn_cast<CallInst>(User)) {
      CallSite CS(CI);
      const Value *calleeV = CS.getCalledValue();
      const Function *callee = dyn_cast<Function>(calleeV->stripPointerCasts());
      if (callee && isVeriCall(*callee)) {
        if (nonBoolCond) {
          FunctionType *FTy = callee->getFunctionType();
          if (!FTy->isVarArg() && FTy->getReturnType()->isVoidTy() &&
              FTy->getNumParams() == 1 && !isBool(FTy->getParamType(0))) {
            verifierCalls.push_back(CI);
            continue;
          }
        } else {
          verifierCalls.push_back(CI);
          continue;
        }
      }
    }
    verifierCalls.clear();
    return false;
  }
  return true;
}

bool AllUsesAreVerifierCalls(Value &V) {
  SmallVector<CallInst *, 4> verifierCalls /*unused*/;
  return AllUsesAreVerifierCalls(V, false, false, verifierCalls);
}

// Return true if all uses are GEPs
bool AllUsesAreGEP(Value &V) {
  for (auto &U : V.uses())
    if (!isa<GetElementPtrInst>(U.getUser()))
      return false;
  return true;
}

bool AllUsesAreIgnoredInst(llvm::Value &V) {
  for (auto &U : V.uses()) {
    if (CallInst *CI = dyn_cast<CallInst>(U.getUser())) {
      if (Function *CalledF = dyn_cast<Function>(CI->getCalledOperand())) {
        if (CalledF->getName().startswith("llvm.dbg.value") ||
            CalledF->getName().startswith("llvm.lifetime")) {
          continue;
        }
      }
    }
    return false;
  }
  return true;
}

} // namespace clam

namespace {
using BasicBlockPtrSet = llvm::SmallPtrSet<const llvm::BasicBlock *, 32>;
} // namespace 

namespace llvm {
template <> class po_iterator_storage<BasicBlockPtrSet, true> {
  BasicBlockPtrSet &Visited;
public:
  po_iterator_storage(BasicBlockPtrSet &VSet) : Visited(VSet) {}
  po_iterator_storage(const po_iterator_storage &S) : Visited(S.Visited) {}
  bool insertEdge(Optional<const BasicBlock *> src, const BasicBlock *dst)
  { return Visited.insert(dst).second; }
  void finishPostorder(const BasicBlock *bb) {}
};
} // namespace llvm

namespace clam {
void revTopoSort(const llvm::Function &F, std::vector<const BasicBlock *> &out) {
  if (F.isDeclaration() || F.empty()) return;
  auto *f = &F;
  BasicBlockPtrSet Visited;
  std::copy(po_ext_begin(f, Visited), po_ext_end(f, Visited), std::back_inserter(out));
}

void topoSort(const llvm::Function &F, std::vector<const BasicBlock*> &out) {
  revTopoSort(F, out);
  std::reverse(out.begin(), out.end());
}

} // namespace clam
